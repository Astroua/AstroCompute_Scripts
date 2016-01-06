############################################################################################################
#CASA Script #2-->Timing script
#Input: Calibrated and Split MS
#Output: Produces a lightcurve with a user specified time bin (plot + data file)
#Note: This script is theoretically compatible with any data that can be imported
#into CASA, but has only been tested on continuum data from the VLA, SMA, and NOEMA
#(import NOEMA data into CASA: http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf).
#############################################################################################################
#Original version written by C. Gough (student of J. Miller-Jones)--> 09/2012
#Last updated by A. Tetarenko--> 10/2015
#############################################################################################################
#Import modules
#
#To ensure all modules used in this script are callable:
#1. download the casa-python executable wrapper package and then you can install any python package to use in CASA
#with the prompt casa-pip --> https://github.com/radio-astro-tools/casa-python (astropy,pyfits,jdcal,photutils,lmfit)
#2. Need uvmultifit -->http://nordic-alma.se/support/software-tools, which needs g++/gcc and gsl libraries
#(http://askubuntu.com/questions/490465/install-gnu-scientific-library-gsl-on-ubuntu-14-04-via-terminal)
#3. Need analysis utilities--> https://casaguides.nrao.edu/index.php?title=Analysis_Utilities
import tempfile
import os
import linecache
import find
from os import path
import numpy as np
import math as m
from jdcal import gcal2jd,jd2gcal
import astropy
from astropy import units as u
from astropy.coordinates import SkyCoord
from datetime import time, timedelta, datetime
import matplotlib.pyplot as pp
from scipy.stats import norm
import re
import sys
import json

from utils import is_power2, load_json

def run_aegean(tables,cellSize):
    '''Loads in and parses data file output from Aegean object detection script (Aegean_ObjDet.py),
    to extract positional information on sources in field

    tables: data file output form Aegean_ObjDet.py
    cellSize: imaging parameter, arcsec/pix

    return: lists of source #, RA, DEC, semi-major axis, semi-minor axis, and position angle
    for all detected sources
    '''
    src_list=[]
    ra_list=[]
    dec_list=[]
    maj_list=[]
    min_list=[]
    pos_list=[]
    cellSize_string=cellSize[0]
    cellSize_list=re.findall('\d+|\D+', cellSize_string)
    cellSize0=float(cellSize_list[0]+cellSize_list[1]+cellSize_list[2])
    with open(tables) as f:
    	lines=f.readlines()
    for i in range(1,len(lines)):
    	lin_split=lines[i].split('\t')
    	src_list.append(lin_split[0])
    	ra_list.append(lin_split[4])#string
    	dec_list.append(lin_split[5])#string
    	maj_list.append(float(lin_split[14])/cellSize0)#pix
    	min_list.append(float(lin_split[16])/cellSize0)#pix
    	pos_list.append(float(lin_split[18]))#deg
    return(src_list,ra_list,dec_list,maj_list,min_list,pos_list)

################################################################################################################
#User Input Section and Setup-->read in from parameters file
#WARNING-->The following variables need to be carefully considered and changed for each new data set
################################################################################################################

#set initial path to where input/output is to be stored
#NOTE: MS's need to be in path_dir/data, all output needs to be in path_dir/data_products,
#both directories need to be created beforehand
# path_dir='/home/ubuntu/'
path_dir = sys.argv[-1]

param_file = sys.argv[-2]

#get input parameters from file
data_params = load_json(param_file)

'''DATA SET PARAMETERS'''
# Target name
target = data_params["target"]
# Date of observation
obsDate = data_params["obsDate"]
# Observation frequency
refFrequency =data_params["refFrequency"]
# Label for casa output directories and files
label = target + '_' + refFrequency + '_' + obsDate + '_'
# Length of time bins (H,M,S); see below if you want manual input (line 288)
intervalSizeH = data_params["intervalSizeH"]
intervalSizeM = data_params["intervalSizeM"]
intervalSizeS = data_params["intervalSizeS"]
# Name of visibliity - should include full path if script is not being run from vis location.
visibility = path_dir+'data/'+ data_params["visibility"]
visibility_uv = path_dir+'data/'+ data_params["visibility"]

'''DIRECTORY AND FILE NAME PARAMETERS'''
# Set path to directory where all output from this script is saved.
outputPath = path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec/'

# dataPath contains the path and filename in which data file will be saved.
# This script can be run on several epochs of data from the same observation without changing this path.
# In this case the data file will be appended each time.
dataPath = path_dir+'data_products/datafile_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec.txt'

#make output directory (within data_products directory)--> done in initial clean, but check if didn't run that
if not os.path.isdir(path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'):
	mkdir_string='sudo mkdir '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
	mkdir_perm1='sudo chown ubuntu '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
	mkdir_perm2='sudo chmod -R 777 '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
	os.system(mkdir_string)
	os.system(mkdir_perm1)
	os.system(mkdir_perm2)

#Object detection file-->output from Aegean_ObjDet.py
tables=outputPath+label+'whole_dataset_objdet_comp.tab'

'''FLAGS'''
#flag to run object detection
runObj=data_params["runObj"]
#Is the data set large enough that you only want to save a cutout?
#If cutout='T' & big_data='T' --> clean,fit, cutout, delete original image each interval
#If cutout='T' & big_data='F' --> clean all, fit all, then delete.
#If cutout='F' & big_data='F' --> clean all full size, fit all full size, no delete
cutout=data_params["cutout"]
big_data=data_params["big_data"]
# Clean can be toggled on/off here (T/F).
runClean =data_params["runClean"]
#optimize cleaning
opt_clean=data_params["opt_clean"]
#do you want to fix parameters in fits from full data set fit? (T of F)
fix_pos=data_params["fix_pos"]
#if fixed parameters do you want to mc sample the fixed parameters (T) or take the best fit (F)?
do_monte=data_params["do_monte"]
#do you want peak (mJy/beam; F) or integrated (mJy; T) flux, or both(B) in lightcurve file?
integ_fit=data_params["integ_fit"]
#do you want to do uv fitting (T or F) as well?
#Source parameters are: x offset (arcsec east), y offset (arcsec north),flux (Jy);
uv_fit=data_params["uv_fit"]
#if fixed parameters do you want to mc sample the fixed parameters (T) or take the best fit (F)?
do_monte_uv=data_params["do_monte_uv"]
#fix position in UV
uv_fix=data_params["uv_fix"]
#If runClean=F set fit_cutout='T' if you have only cutout images but want to refit without cleaning again,
fit_cutout=data_params["fit_cutout"]
#define start and end time yourself
def_times=data_params["def_times"]

'''CLEAN PARAMETERS'''
# The clean command (line 505) should be inspected closely to ensure all arguments are appropriate before
# running script on a new data set.
# The following arguments will be passed to casa's clean, imfit or imstat functions:
imageSize = [data_params["imageSize"]* 2
if opt_clean=='T':
	if is_power2(imageSize[0])==False:
		print 'Clean will run faster if image size is 2^n'
		imageSize=[int(pow(2, m.ceil(np.log(imageSize[0])/np.log(2)))),int(pow(2, m.ceil(np.log(imageSize[0])/np.log(2))))]
        	print 'imagesize is now set to ', imageSize
    	print 'imagesize remains at ',imageSize, 'due to user request'
numberIters = data_params["numberIters"]
cellSize = [data_params["cellSize"] * 2
taylorTerms = data_params["taylorTerms"]
myStokes = data_params["myStokes"]
thre=data_params["thre"]
spw_choice=data_params["spw_choice"]
# If an outlier file is to be used in the clean, set outlierFile to the filename (path inluded). myThreshold will
# also need to be set if outlier file is to be used.
# If not, set outlierFile to ''.
outlierFile = data_params["outlierFile"]


'''OBJECT DETECTION AND SELECTION PARAMETERS'''
if runObj=='T':
    #object detection with Aegean algorithm--> Need to run initial_clean.py in CASA, and Aegean_ObjDet.py outside
    #CASA first
    src_l,ra_l,dec_l,maj_l,min_l,pos_l=run_aegean(tables,cellSize)
    ind=data_params["ind"]
    # target position-->Take bounding ellipse from Aegean and convert to minimum bounding box in pixels for
    #use with rest of script
    tar_pos=au.findRADec(outputPath+label+'whole_dataset.image',ra_l[int(ind)-1]+' '+dec_l[int(ind)-1])
    bbox_halfwidth=np.sqrt((min_l[int(ind)-1]*np.cos(pos_l[int(ind)-1]))**2+(min_l[int(ind)-1]*np.sin(pos_l[int(ind)-1]))**2)+3
    bbox_halfheight=np.sqrt((maj_l[int(ind)-1]*np.cos(pos_l[int(ind)-1]+(np.pi/2.)))**2+(maj_l[int(ind)-1]*np.sin(pos_l[int(ind)-1]+(np.pi/2.)))**2)+3
    targetBox = str(tar_pos[0]-bbox_halfwidth)+','+ str(tar_pos[1]-bbox_halfheight)+','+str(tar_pos[0]+bbox_halfwidth)+','+ str(tar_pos[1]+bbox_halfheight)
elif runObj=='F':
#input target box in pixels if not running object detection
    targetBox =data_params["targetBox"] # #'2982,2937,2997,2947'
#mask for clean based on target box
maskPath = 'box [['+targetBox.split(',')[0]+'pix,'+targetBox.split(',')[1]+'pix],['+targetBox.split(',')[2]+'pix,'+targetBox.split(',')[3]+'pix]]'#path_dir+'data/v404_jun22B_K21_clean_psc1.mask'

'''IMAGE PRODUCT PARAMETERS: CUTOUT AND RMS/ERROR IN IMAGE PLANE HANDLING'''
#define rms boxes for realistic error calculation
#pix_shift_cutout is how many pixels past target bx you want in cutout images rms calculation
pix_shift_cutout=data_params["pix_shift_cutout"]
annulus_rad_inner=data_params["annulus_rad_inner"]
annulus_rad_outer=data_params["annulus_rad_outer"]
cut_reg=str(float(targetBox.split(',')[0])-pix_shift_cutout)+','+str(float(targetBox.split(',')[1])-pix_shift_cutout)+','+str(float(targetBox.split(',')[2])+pix_shift_cutout)+','+str(float(targetBox.split(',')[3])+pix_shift_cutout)#'2962,2917,3017,2967'
x_sizel=float(targetBox.split(',')[0])
x_sizeu=float(targetBox.split(',')[2])
y_sizel=float(targetBox.split(',')[1])
y_sizeu=float(targetBox.split(',')[3])
#local rms around target box but inside cutout region (not as accurate as annulus)
#rmsbox1=str(x_sizel-(3./4.)*pix_shift_cutout)+','+str(y_sizel-(3./4.)*pix_shift_cutout)+','+str(x_sizel-(1./4.)*pix_shift_cutout)+','+str(y_sizel-(1./4.)*pix_shift_cutout)
#rmsbox2=str(x_sizeu+(1./4.)*pix_shift_cutout)+','+str(y_sizel-(3./4.)*pix_shift_cutout)+','+str(x_sizeu+(3./4.)*pix_shift_cutout)+','+str(y_sizel-(1./4.)*pix_shift_cutout)
#rmsbox3=str(x_sizel-(1./4.)*pix_shift_cutout)+','+str(y_sizeu+(1./4.)*pix_shift_cutout)+','+str(x_sizeu+(1./4.)*pix_shift_cutout)+','+str(y_sizeu+(3./4.)*pix_shift_cutout)
#annulus region parameters
cen_annulus='['+str(((x_sizeu-x_sizel)/2.)+x_sizel)+','+str(((y_sizeu-y_sizel)/2.)+y_sizel)+']'
cen_radius='['+str(annulus_rad_inner)+'pix,'+str(annulus_rad_outer)+'pix]'
#what unit do you want light curve
lc_scale_unit=data_params["lc_scale_unit"]
if lc_scale_unit=='m':
	lc_scale_factor=1.0e03
elif lc_scale_unit=='u':
	lc_scale_factor=1.0e06
elif lc_scale_unit=='':
	lc_scale_factor=1.0
#what time scale do you want in light curve
lc_scale_time=data_params["lc_scale_time"]
#if remainder of obstime/interval > this # then it is included in time intervals (in minutes)
rem_int=data_params["rem_int"]

'''REFITTTING WITHOUT CLEANING PARAMETERS'''
#make sure rms boxes are set to local.
if fit_cutout=='T':
	targetBox = str(float(targetBox.split(',')[0])-float(cut_reg.split(',')[0]))+','+str(float(targetBox.split(',')[1])-float(cut_reg.split(',')[1]))+','+str(float(targetBox.split(',')[2])-float(cut_reg.split(',')[0]))+','+str(float(targetBox.split(',')[3])-float(cut_reg.split(',')[1]))
	cx_sizel=float(targetBox.split(',')[0])
	cx_sizeu=float(targetBox.split(',')[2])
	cy_sizel=float(targetBox.split(',')[1])
	cy_sizeu=float(targetBox.split(',')[3])
	#rmsbox1=str(cx_sizel-(3./4.)*pix_shift_cutout)+','+str(cy_sizel-(3./4.)*pix_shift_cutout)+','+str(cx_sizel-(1./4.)*pix_shift_cutout)+','+str(cy_sizel-(1./4.)*pix_shift_cutout)
	#rmsbox2=str(cx_sizeu+(1./4.)*pix_shift_cutout)+','+str(cy_sizel-(3./4.)*pix_shift_cutout)+','+str(cx_sizeu+(3./4.)*pix_shift_cutout)+','+str(cy_sizel-(1./4.)*pix_shift_cutout)
	#rmsbox3=str(cx_sizel-(1./4.)*pix_shift_cutout)+','+str(cy_sizeu+(1./4.)*pix_shift_cutout)+','+str(cx_sizeu+(1./4.)*pix_shift_cutout)+','+str(cy_sizeu+(3./4.)*pix_shift_cutout)
	cen_annulus='['+str(((cx_sizeu-cx_sizel)/2.)+cx_sizel)+','+str(((cy_sizeu-cy_sizel)/2.)+cy_sizel)+']'
	cen_radius='[10pix,20pix]'

'''IMAGE FITTING PARAMETERS'''
if do_monte == 'T':
#add appropriate label to final lightcurve file name
	lab='_mc_'
else:
	lab='_bestfit_'
# number of MC simulations if MC position sampling chosen
nsim=100
#if fixing parameters, what parameters do you want to fix in fits? f= peak flux, x=peak x pos, y=peak y pos, a=major axis (arcsec), b=minor axis (arcsec), p=position angle (deg)
#a, b, p convolved with beam values not deconvolved values!!!
#format is [value,error] from fit of full data set
if fix_pos=='T':
    par_fix=data_params["par_fix"]
    if runObj=='F':
    	print 'Cleaning Full Data Set-->'
    	clean(vis=visibility, imagename=outputPath+label+'whole_dataset', field='', mode='mfs', imsize=imageSize, cell=cellSize, weighting='natural',spw=spw_choice, nterms=taylorTerms, niter=numberIters, gain=0.1, threshold=thre, interactive=F)
    print 'Fitting full data set in Image Plane-->'
    full_fit=imfit(imagename=outputPath+label+'whole_dataset.image',box=targetBox,logfile=outputPath+label+'whole_dataset.txt')
    imfitFilefull=open(outputPath+label+'whole_dataset.txt','r')
    for line in imfitFilefull:
        if ('--- ra:' in line) & ('pixels' in line):
            ra_string=line
        if ('--- dec:' in line) & ('pixels' in line):
            dec_string=line
    pmra=ra_string.find('+/-')
    pmdec=dec_string.find('+/-')
    pixra=ra_string.find('pixels')
    pixdec=dec_string.find('pixels')

    peak_x=[float(ra_string[16:pmra]),float(ra_string[29:pixra])]#[2988.63,0.02]
    peak_y=[float(dec_string[16:pmdec-1]),float(dec_string[29:pixdec-1])]#[2942.09,0.01]
    if fit_cutout=='T':
        peak_x=[peak_x[0]-float(cut_reg.split(',')[0]),peak_x[1]]
        peak_y=[peak_y[0]-float(cut_reg.split(',')[1]),peak_y[1]]
    b_maj=[full_fit['results']['component0']['shape']['majoraxis']['value'],full_fit['results']['component0']['shape']['majoraxiserror']['value']]#[0.154,0.001]
    b_min=[full_fit['results']['component0']['shape']['minoraxis']['value'],full_fit['results']['component0']['shape']['minoraxiserror']['value']]#[0.099,0.0005]
    pos_ang=[full_fit['results']['component0']['shape']['positionangle']['value'],full_fit['results']['component0']['shape']['positionangleerror']['value']]#[67.41,0.42]

'''UV FITTING PARAMETERS'''
#only point sources right now
comp_uv='delta'
#Not tested on anything other than I; for VLA 'I', for SMA 'LL' (from listobs)
stokes_param=data_params["stokes_param"]
if uv_fit=='T':
    print 'Fitting full data set in UV Plane-->'
    fitfulluv=uvm.uvmultifit(vis=visibility_uv, spw=spw_choice, column = "data", uniform=False, write_model=False, model=[comp_uv],stokes = stokes_param, var=['p[0],p[1],p[2]'],outfile = outputPath+label+'whole_dataset_uv.txt', OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
    if uv_fix=='F':
    	uv_var='p[0],p[1],p[2]'#'2.8194e-02,8.5502e-03,p[2]'
    elif uv_fix=='T':
    	uv_var=str(fitfulluv.result['Parameters'][0])+','+str(fitfulluv.result['Parameters'][1])+',p[2]'#'2.8194e-02,8.5502e-03,p[2]'
    src_uv_init=str(fitfulluv.result['Parameters'][0])+','+str(fitfulluv.result['Parameters'][1])+','+str(fitfulluv.result['Parameters'][2])#[2.8194e-02,8.5502e-03 , 1.3508e-01]
    src_uv_err=str(fitfulluv.result['Uncertainties'][0])+','+str(fitfulluv.result['Uncertainties'][1])+','+str(fitfulluv.result['Uncertainties'][2])#[4.7722e-05 , 3.7205e-05, 1.1192e-04]

############################################################################################################
#End of User Input Section and Setup
############################################################################################################

############################################################################################################
#Calculate time bins
############################################################################################################

# Run listobs and store as a text file.
listobs(vis=visibility, listfile=outputPath + label + 'listobs.text')

# Get time range of observation from listobs and generate timeIntervals vector.
# Extract start and end times and start date from listobs. The first few lines of listobs output
# are identical between listobs executions. Therefore we can be confident that the time info will
# always be located at the same position.

listobsLine7 = linecache.getline(outputPath + label + 'listobs.text', 7)

startTimeH = int(listobsLine7[31:33])
startTimeM = int(listobsLine7[34:36])
startTimeS= int(listobsLine7[37:39])
endTimeH = int(listobsLine7[61:63])
endTimeM = int(listobsLine7[64:66])
endTimeS = int(listobsLine7[67:69])
startDate = listobsLine7[19:30]
endDate = listobsLine7[49:60]

# The data will be plotted against MJD time, so we need startDate in this format also.
# Convert month to integer.
year = startDate[7:11]
month = startDate[3:6]
day = startDate[0:2]
monthInt = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}[month]
startDateMJD = gcal2jd(year,monthInt,day)
yeare = endDate[7:11]
monthe = endDate[3:6]
daye = endDate[0:2]
monthInte = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}[monthe]
endDateMJD = gcal2jd(yeare,monthInte,daye)

####################################################################################################
# Optional: Define start and end times manually (i.e., if only want subset of observation )
#
if def_times=='T':
	startTimeH = data_params["startTimeH"] #input('Enter start time hour >| ')
	startTimeM = data_params["startTimeM"] #input('Enter start time minute >| ')
	startTimeS = data_params["startTimeS"] #input('Enter start time seconds >| ')
	endTimeH = data_params["endTimeH"] #input('Enter finish time hour >| ')
	endTimeM = data_params["endTimeM"] #input('Enter finish time minute >| ')
	endTimeS = data_params["endTimeS"] #input('Enter start time seconds >| ')
####################################################################################################

# We require a list of time interval strings to pass to the clean function.
# In order to generate this list we need to perform some basic arithmetic on time values.
# The datetime module achieves this using timedelta objects.

# First the start and end times are converted to datetime.time objects.
startTime = time(startTimeH,startTimeM,startTimeS)
endTime = time(endTimeH,endTimeM,endTimeS)
#the following is used to enable input of multi-day observations, essentially both dates and
#times are recorded
testTime = time(23,59,59)
test2Time = time(00,00,01)
#
# Next we convert these to datetime.timedelta objects.
startTimeDelta = timedelta(hours=startTime.hour, minutes=startTime.minute,seconds=startTime.second)
endTimeDelta = timedelta(hours=endTime.hour, minutes=endTime.minute,seconds=endTime.second)
temp = timedelta(hours=testTime.hour, minutes=testTime.minute,seconds=testTime.second)
temp2 = timedelta(hours=test2Time.hour, minutes=test2Time.minute,seconds=test2Time.second)
#
# Use timedelta to calculate duration.
if float(startDateMJD[1]) == float(endDateMJD[1]):
	observationDurationDelta = endTimeDelta-startTimeDelta
else:
	observationDurationDelta = endTimeDelta+(temp-startTimeDelta+temp2)
# The result is converted back to a datetime.time object.
observationDuration = (datetime.min+observationDurationDelta).time()
#
# Duration is printed in iso format as a string for the user.
print '\nTotal observation time is ' + str(time.isoformat(observationDuration))
#
# User is prompted for the desired interval length.
#intervalSizeH = input('Enter length of interval (hours) >| ')
#intervalSizeM = input('Enter length of interval (minutes) >| ')
#intervalSizeS = input('Enter length of interval (seconds) >| ')
#
# The interval length is converted to a datetime.timedelta object
intervalSize = time(intervalSizeH, intervalSizeM, intervalSizeS)
intervalSizeDelta = timedelta(hours=intervalSizeH, minutes=intervalSizeM, seconds=intervalSizeS)
#
# Calculate number of intervals. This is done by converting the duration and interval size, which are both datetime.time objects, into an integer number of minutes. Integer division is then used.
durationSeconds = (observationDuration.hour)*3600 + (observationDuration.minute)*60+(observationDuration.second)
intervalSeconds = (intervalSize.hour)*3600 + (intervalSize.minute)*60 + (intervalSize.second)
numIntervals = durationSeconds / intervalSeconds
print numIntervals

# The remaining observation time is calculated. This will be used if it is > rem_int.
remainder = durationSeconds % intervalSeconds

# The time intervals are generated by a for loop as a list of strings.
# These are now ready to be passed to CLEAN. A list of MJD times is generated simultaneously,
# which will be used to label the light curve. The MJD times refer to the end time of each interval.

# The lists are initialised.
timeIntervals = range(numIntervals)
mjdTimes = range(numIntervals)

for element in range(numIntervals):
    # Start and end times must be re-initialised at start of each iteration.
    startTimeDelta = timedelta(hours=startTime.hour, minutes=startTime.minute,seconds=startTime.second)
    endTimeDelta = timedelta(hours=endTime.hour, minutes=endTime.minute,seconds=endTime.second)
    # This for loop adds x+1 intervals to the start time, where x is the number of iterations
    # perfomed by the outer loop so far.
    for x in range(element):
        startTimeDelta = startTimeDelta + intervalSizeDelta
    # The start time of the interval should now be correct, so the end time is found by adding
    # one more interval length to the start time.
    endTimeDelta = startTimeDelta + intervalSizeDelta
    # As mentioned above, the MJD times list lists only the end time of each interval (in MJD format).
    # The next element in the list is generated by adding the fractional number of days since the start
    # of the observation to the MJD start time.
    if float(endTimeDelta.days) == 0:
    	mjdTimes[element] = startDateMJD[1] + long((datetime.min+endTimeDelta).time().hour)/24.0 + long((datetime.min+endTimeDelta).time().minute)/60.0/24.0 + long((datetime.min+endTimeDelta).time().second)/60.0/60.0/24.0
    	date_conv=jd2gcal(startDateMJD[0],startDateMJD[1])
    	month_0={1: '01', 2: '02', 3: '03', 4: '04', 5: '05', 6: '06', 7: '07', 8: '08', 9: '09', 10: '10', 11: '11', 12: '12'}[date_conv[1]]

    else:
    	mjdTimes[element] = endDateMJD[1] + long((datetime.min+endTimeDelta).time().hour)/24.0 + long((datetime.min+endTimeDelta).time().minute)/60.0/24.0 + long((datetime.min+endTimeDelta).time().second)/60.0/60.0/24.0
    	date_conv2=jd2gcal(endDateMJD[0],endDateMJD[1])
    	month_2={1: '01', 2: '02', 3: '03', 4: '04', 5: '05', 6: '06', 7: '07', 8: '08', 9: '09', 10: '10', 11: '11', 12: '12'}[date_conv2[1]]

    # The start and end times of the interval are converted to strings in a format ready to be
    # passed to clean, ie. 'hh:mm:ss~hh:mm:ss'.
    if float(endTimeDelta.days) == 0 and float(startTimeDelta.days)==0:
    	timeIntervals[element] = str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+startTimeDelta).time())) + '~' +str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+endTimeDelta).time()))
    elif float(endTimeDelta.days) != 0 and float(startTimeDelta.days)==0:
    	timeIntervals[element] = str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+startTimeDelta).time())) + '~' +str(date_conv2[0]) +'/'+ str(month_2) +'/'+str(str(date_conv2[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+endTimeDelta).time()))
    elif float(endTimeDelta.days) != 0 and float(startTimeDelta.days)!=0:
    	timeIntervals[element] = str(date_conv2[0]) +'/'+ str(month_2) +'/'+str(str(date_conv2[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+startTimeDelta).time())) + '~' +str(date_conv2[0]) +'/'+ str(month_2) +'/'+str(str(date_conv2[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+endTimeDelta).time()))
    else:
    	print 'Something is wrong, please review input'


# If the remainder of observation time is greater than this # of minutes it should be used, and is appended
# to timeIntervals and mjdTimes.
if remainder >= rem_int:
    if float(startDateMJD[1]) != float(endDateMJD[1]):
    	timeIntervals = timeIntervals + [str(date_conv2[0])+'/'+str(month_2)+'/'+str(str(date_conv2[2]).zfill(2)) +'/'+str(time.isoformat((datetime.min+endTimeDelta).time()))+'~'+str(date_conv2[0])+'/'+str(month_2)+'/'+str(str(date_conv2[2]).zfill(2)) +'/'+str(time.isoformat(endTime))]
    	mjdTimes = mjdTimes + [endDateMJD[1] + long(endTime.hour)/24.0 + long(endTime.minute)/24.0/60.0+ long(endTime.second)/24.0/60.0/60.0]
    else:
    	timeIntervals = timeIntervals + [str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +'/' +str(time.isoformat((datetime.min+endTimeDelta).time()))+'~'+str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +'/' +str(time.isoformat(endTime))]
    	mjdTimes = mjdTimes + [startDateMJD[1] + long(endTime.hour)/24.0 + long(endTime.minute)/24.0/60.0+ long(endTime.second)/24.0/60.0/60.0]

# The results are printed for the user.
print '\nThe observation will be divided into the following intervals: '
print timeIntervals
print '\nmjdTimes'
print mjdTimes

############################################################################################
# Clean each chunk individually and fit.
############################################################################################
# Imfit will be run on each image with a given set of estimates. These estimates need to be
# given as a text file with a particular layout.
# The estimates can be gathered from the output of running imhead() on each image.
# For each image the for loop runs imhead(), creates a temporary estimates text file from its
#output, and then runs imfit() using these estimates. The resulting imfits are saved as
#text files.
############################################################################################

#initialize lists
#result_box1rms=[]
#result_box2rms=[]
#result_box3rms=[]
fluxError_real=[]
uv_fitval=[]
uv_fiterr=[]
#make copy of mjdTimes and timeIntervals for uv fitting
mjdTimes_uv=mjdTimes
timeIntervals_uv=timeIntervals

print 'Clean is starting-->'
if runClean:
	for interval, time, interval_uv, time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
		print 'cleaning interval:', interval
        	if outlierFile == '':
            		intervalString=interval.replace(':', '.').replace('/','_')
            		clean(vis=visibility, imagename=outputPath+label+intervalString, mask=maskPath, selectdata=T, timerange=interval, field='', mode='mfs', imsize=imageSize, cell=cellSize, weighting='natural', usescratch=T,spw=spw_choice, nterms=taylorTerms, niter=numberIters, gain=0.1, threshold=thre, interactive=F)

			if big_data == 'T':

				#

		    		intervalString=interval.replace(':', '.').replace('/','_')
		    		print intervalString,

		   		 # Depending on the nuber of terms used for the Taylor expansion in clean, the image file name will end in either
		    		# .image or .image.tt0. This is handled with the following if statement.
		    		if outlierFile != '':
					imSuffix = '_0.image'
		    		elif taylorTerms == 1:
					imSuffix = '.image'
		    		else:
					imSuffix = '.image.tt0'

		    		# For some intervals the CLEAN may have failed, in which case the image file for that interval will not exist.
		   		 # An if statement is used to pass over these intervals, avoiding runtime errors.
		    		if os.path.exists(outputPath+label+intervalString+imSuffix):
					# In most cases imhead returns a dictionary with the value located by the key 'value'
					beamMajor = imhead(imagename=outputPath+label+intervalString+imSuffix,mode='get',hdkey='beammajor')
					beamMajor = str(beamMajor['value'])+'arcsec'
					beamMinor = imhead(imagename=outputPath+label+intervalString+imSuffix,mode='get',hdkey='beamminor')
					beamMinor = str(beamMinor['value'])+'arcsec'
					beamPA = imhead(imagename=outputPath+label+intervalString+imSuffix,mode='get',hdkey='beampa')
					beamPA = str(beamPA['value'])+'deg'
					imstatOut = imstat(imagename=outputPath+label+intervalString+imSuffix,box=targetBox)
					peak = imstatOut['max']
					peak = str(peak[0])
					peakPosValue = imstatOut['maxpos']
					# In this case 'value' references an array. The x and y coordinates can be retieved by indexing.
					peakPosX = str(peakPosValue[0])
					peakPosY = str(peakPosValue[1])
					# Save parameters in a temp file which will be passed to imfit. The correct layout is:
					# peak intensity, peak x-pixel value, peak y-pixel value, major axis, minor axis, position angle, fixed
					# the fixed parameter can contain any of the following:
					#'f' (peak intensity), 'x' (peak x position), 'y' (peak y position), 'a' (major axis), 'b' (minor axis), 'p' (position angle)
		#        tempFile = tempfile.NamedTemporaryFile()
					tempFile = open('tempfile.txt','w')
					mystring = str(peak+', '+peakPosX+', '+peakPosY+', '+beamMajor+', '+beamMinor+', '+beamPA+',abp')
					tempFile.write(mystring)
		    			# Run imfit using the above beam parameters.
					tempFile.close()
					print timeIntervals.index(interval),':',outputPath+label+intervalString+imSuffix
					imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text', estimates=tempFile.name, append=F, overwrite = T)
		    			if fix_pos == 'T':
						if do_monte == 'T':
							samp_px=np.random.normal(0,1,nsim)
							samp_py=np.random.normal(0,1,nsim)
							samp_bmaj=np.random.normal(0,1,nsim)
							samp_bmin=np.random.normal(0,1,nsim)
							samp_pos=np.random.normal(0,1,nsim)
							for i in range(0,len(samp_px)):
								peak_x1=(samp_px[i]*peak_x[1])+peak_x[0]
								peak_y1=(samp_py[i]*peak_y[1])+peak_y[0]
								b_maj1=(samp_bmaj[i]*b_maj[1])+b_maj[0]
								b_min1=(samp_bmin[i]*b_min[1])+b_min[0]
								pos_ang1=(samp_pos[i]*pos_ang[1])+pos_ang[0]
								tempFile2 = open('tempfile2.txt','w')
								mystring2 = str(peak+', '+str(peak_x1)+', '+str(peak_y1)+', '+str(b_maj1)+'arcsec, '+str(b_min1)+'arcsec, '+str(pos_ang1)+'deg, '+par_fix)
								tempFile2.write(mystring2)
								tempFile2.close()
								imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text', estimates=tempFile2.name,append=F, overwrite = T)
						elif do_monte =='F':
		    					tempFile2 = open('tempfile2.txt','w')
							mystring2 = str(peak+', '+str(peak_x[0])+', '+str(peak_y[0])+', '+str(b_maj[0])+'arcsec, '+str(b_min[0])+'arcsec, '+str(pos_ang[0])+'deg,'+par_fix)
							tempFile2.write(mystring2)
		    					tempFile2.close()
							imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text', estimates=tempFile2.name,append=F, overwrite = T)
						else:
							print 'Please specify whether you wish to perform a Monte Carlo fit (T) or not(F)'
		    			else:
		    				imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text',append=F, overwrite = T)

					result_box1=imstat(imagename=outputPath+label+intervalString+imSuffix,region='annulus['+cen_annulus+','+cen_radius+']')
					#result_box2=imstat(imagename=outputPath+label+intervalString+imSuffix,box=rmsbox2)
					#result_box3=imstat(imagename=outputPath+label+intervalString+imSuffix,box=rmsbox3)
					#result_box1rms.append(result_box1['rms'][0])
					#result_box2rms.append(result_box2['rms'][0])
					#result_box3rms.append(result_box3['rms'][0])
					#rms_array=np.array([result_box1['rms'],result_box2['rms'],result_box3['rms']])
					fluxError_real.append(result_box1['rms'][0])
		    		else:
					print '\nCLEAN failed on interval ' + interval + '.'
					# The corresponding time intervals must be removed from timeIntervals to avoid runtime errors further on.
					timeIntervals.remove(interval)
					mjdTimes.remove(time)
					#timeIntervals_uv.remove(interval_uv)
					#mjdTimes_uv.remove(time_uv)
		    		if uv_fit =='T' and os.path.exists(outputPath+label+intervalString+imSuffix):
					if do_monte_uv == 'T':
						monte_uvval=[]
						samp_x_uv=np.random.normal(0,1,nsim)
						samp_y_uv=np.random.normal(0,1,nsim)
						for i in range(0,len(samp_x_uv)):
							x_uv1=(samp_x_uv[i]*src_uv_err[1])+src_uv[1]
							y_uv1=(samp_y_uv[i]*src_uv_err[2])+src_uv[2]
							if np.where(np.array(timeIntervals)==interval)[0][0]==0:
								fit=uvm.uvmultifit(vis=visibility_uv, MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv],spw=spw_choice, column = "data", uniform=False, write_model=False, model=[comp_uv],stokes = stokes_param, var=[str(x_uv1)+','+str(y_uv1)+','+'p[2]'], p_ini=[x_uv1,y_uv1,src_uv_init[2]],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'_'+str(i)+'.txt', OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
							else:
								fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
								fit.var=[str(x_uv1)+','+str(y_uv1)+','+'p[2]']
								fit.p_ini=[x_uv1,y_uv1,src_uv_init[2]]
								fit.fit(redo_fixed=False,reinit_model=False)
							monte_uvval.append(fit.result['Parameters'][2])
						uv_fitval.append(monte_uvval)

					elif do_monte_uv =='F':
						if np.where(np.array(timeIntervals)==interval)[0][0]==0:
							fit=uvm.uvmultifit(vis=visibility_uv, MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv],spw=spw_choice, column = "data", uniform=False, write_model=False, model=[comp_uv],stokes = stokes_param, var=[uv_var], p_ini=[src_uv_init[0],src_uv_init[1],src_uv_init[2]],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.txt', OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
						else:
							fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
							fit.fit(redo_fixed=False,reinit_model=False)
						uv_fitval.append(fit.result['Parameters'][2])
						uv_fiterr.append(fit.result['Uncertainties'][2])
					else:
						print 'Please specify whether you wish to perform a Monte Carlo fit o uv, (T) or not(F)'
	    			if cutout== 'T':#and os.path.exists(outputPath+label+intervalString+imSuffix):box=rmsbox1
	    				immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',box=cut_reg,outfile=outputPath+label+intervalString+'_temp.image')
					immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',region='annulus['+cen_annulus+','+cen_radius+']',outfile=outputPath+label+intervalString+'_rms.image')
					#immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',box=rmsbox2,outfile=outputPath+label+intervalString+'_rms2.image')
					#immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',box=rmsbox3,outfile=outputPath+label+intervalString+'_rms3.image')

					comm_and1='rm -rf '+outputPath+label+intervalString+'.*'
					comm_and2='mv '+outputPath+label+intervalString+'_temp.image '+outputPath+label+intervalString+'.image'
					os.system(comm_and1)
					os.system(comm_and2)


if big_data=='F' or runClean==F:
	for interval, time, interval_uv,time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
	    # Get beam parameters using imhead.

	    intervalString=interval.replace(':', '.').replace('/','_')
	    print intervalString,

	    # Depending on the nuber of terms used for the Taylor expansion in clean, the image file name will end in either
	    # .image or .image.tt0. This is handled with the following if statement.
	    if outlierFile != '':
		imSuffix = '_0.image'
	    elif taylorTerms == 1:
		imSuffix = '.image'
	    else:
		imSuffix = '.image.tt0'

	    # For some intervals the CLEAN may have failed, in which case the image file for that interval will not exist.
	    # An if statement is used to pass over these intervals, avoiding runtime errors.
	    if os.path.exists(outputPath+label+intervalString+imSuffix):
		# In most cases imhead returns a dictionary with the value located by the key 'value'
		beamMajor = imhead(imagename=outputPath+label+intervalString+imSuffix,mode='get',hdkey='beammajor')
		beamMajor = str(beamMajor['value'])+'arcsec'
		beamMinor = imhead(imagename=outputPath+label+intervalString+imSuffix,mode='get',hdkey='beamminor')
		beamMinor = str(beamMinor['value'])+'arcsec'
		beamPA = imhead(imagename=outputPath+label+intervalString+imSuffix,mode='get',hdkey='beampa')
		beamPA = str(beamPA['value'])+'deg'
		imstatOut = imstat(imagename=outputPath+label+intervalString+imSuffix,box=targetBox)
		peak = imstatOut['max']
		peak = str(peak[0])
		peakPosValue = imstatOut['maxpos']
		# In this case 'value' references an array. The x and y coordinates can be retieved by indexing.
		peakPosX = str(peakPosValue[0])
		peakPosY = str(peakPosValue[1])
		# Save parameters in a temp file which will be passed to imfit. The correct layout is:
		# peak intensity, peak x-pixel value, peak y-pixel value, major axis, minor axis, position angle, fixed
		# the fixed parameter can contain any of the following:
		#'f' (peak intensity), 'x' (peak x position), 'y' (peak y position), 'a' (major axis), 'b' (minor axis), 'p' (position angle)
	#        tempFile = tempfile.NamedTemporaryFile()
		tempFile = open('tempfile.txt','w')
		mystring = str(peak+', '+peakPosX+', '+peakPosY+', '+beamMajor+', '+beamMinor+', '+beamPA+',abp')
		tempFile.write(mystring)
	    	# Run imfit using the above beam parameters.
		tempFile.close()
		print timeIntervals.index(interval),':',outputPath+label+intervalString+imSuffix
		imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text', estimates=tempFile.name, append=F, overwrite = T)
	    	if fix_pos == 'T':
			if do_monte == 'T':
				samp_px=np.random.normal(0,1,nsim)
				samp_py=np.random.normal(0,1,nsim)
				samp_bmaj=np.random.normal(0,1,nsim)
				samp_bmin=np.random.normal(0,1,nsim)
				samp_pos=np.random.normal(0,1,nsim)
				for i in range(0,len(samp_px)):
					peak_x1=(samp_px[i]*peak_x[1])+peak_x[0]
					peak_y1=(samp_py[i]*peak_y[1])+peak_y[0]
					b_maj1=(samp_bmaj[i]*b_maj[1])+b_maj[0]
					b_min1=(samp_bmin[i]*b_min[1])+b_min[0]
					pos_ang1=(samp_pos[i]*pos_ang[1])+pos_ang[0]
					tempFile2 = open('tempfile2.txt','w')
					mystring2 = str(peak+', '+str(peak_x1)+', '+str(peak_y1)+', '+str(b_maj1)+'arcsec, '+str(b_min1)+'arcsec, '+str(pos_ang1)+'deg, '+par_fix)
					tempFile2.write(mystring2)
					tempFile2.close()
					imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text', estimates=tempFile2.name,append=F, overwrite = T)
			elif do_monte =='F':
	    			tempFile2 = open('tempfile2.txt','w')
				mystring2 = str(peak+', '+str(peak_x[0])+', '+str(peak_y[0])+', '+str(b_maj[0])+'arcsec, '+str(b_min[0])+'arcsec, '+str(pos_ang[0])+'deg,'+par_fix)
				tempFile2.write(mystring2)
	    			tempFile2.close()
				imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text', estimates=tempFile2.name,append=F, overwrite = T)
			else:
				print 'Please specify whether you wish to perform a Monte Carlo fit (T) or not(F)'
	    	else:
	    		imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox, logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text',append=F, overwrite = T)
		if fit_cutout=='T':
			result_box1=imstat(imagename=outputPath+label+intervalString+'_rms'+imSuffix,region='annulus['+cen_annulus+','+cen_radius+']')
			#result_box2=imstat(imagename=outputPath+label+intervalString+'_rms2'+imSuffix,box=rmsbox2)
			#result_box3=imstat(imagename=outputPath+label+intervalString+'_rms3'+imSuffix,box=rmsbox3)
		else:
			result_box1=imstat(imagename=outputPath+label+intervalString+imSuffix,box=region='annulus['+cen_annulus+','+cen_radius+']')
			#result_box2=imstat(imagename=outputPath+label+intervalString+imSuffix,box=rmsbox2)
			#result_box3=imstat(imagename=outputPath+label+intervalString+imSuffix,box=rmsbox3)
		#result_box1rms.append(result_box1['rms'][0])
		#result_box2rms.append(result_box2['rms'][0])
		#result_box3rms.append(result_box3['rms'][0])
		#rms_array=np.array([result_box1['rms'],result_box2['rms'],result_box3['rms']])
		fluxError_real.append(result_box1['rms'][0])
		#comm_and1='rm -rf '+outputPath+label+intervalString+'.*'
		#comm_and2='mv temp_'+outputPath+label+intervalString+'.image '+outputPath+label+intervalString+'.image'
		#os.system(comm_and1)
		#os.system(comm_and2)
	    else:
		print '\nCLEAN failed on interval ' + interval + '.'
		# The corresponding time intervals must be removed from timeIntervals to avoid runtime errors further on.
		timeIntervals.remove(interval)
		mjdTimes.remove(time)
		#timeIntervals_uv.remove(interval)
		#mjdTimes_uv.remove(time)
	    if uv_fit =='T' and os.path.exists(outputPath+label+intervalString+imSuffix):
		if do_monte_uv == 'T':
			monte_uvval=[]
			samp_x_uv=np.random.normal(0,1,nsim)
			samp_y_uv=np.random.normal(0,1,nsim)
			for i in range(0,len(samp_x_uv)):
				x_uv1=(samp_x_uv[i]*src_uv_err[1])+src_uv[1]
				y_uv1=(samp_y_uv[i]*src_uv_err[2])+src_uv[2]
				if np.where(np.array(timeIntervals)==interval)[0][0]==0:
					fit=uvm.uvmultifit(vis=visibility_uv, MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv],spw=spw_choice, column = "data", uniform=False, write_model=False, model=[comp_uv],stokes = stokes_param, var=[str(x_uv1)+','+str(y_uv1)+','+'p[2]'], p_ini=[x_uv2,y_uv1,src_uv_init[2]],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'_'+str(i)+'.txt', OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
				else:
					fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
					fit.var=[str(x_uv1)+','+str(y_uv1)+','+'p[2]']
					fit.p_ini=[x_uv2,y_uv1,src_uv_init[2]]
					fit.fit(redo_fixed=False,reinit_model=False)
				monte_uvval.append(fit.result['Parameters'][2])
			uv_fitval.append(monte_uvval)
		elif do_monte_uv =='F':
			if np.where(np.array(timeIntervals)==interval)[0][0]==0:
				fit=uvm.uvmultifit(vis=visibility_uv,MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv], spw=spw_choice, column = "data", uniform=False, write_model=False, model=[comp_uv],stokes = stokes_param, var=[uv_var], p_ini=[src_uv_init[0],src_uv_init[1],src_uv_init[2]],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.txt', OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
			else:
				fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
				fit.fit(redo_fixed=False,reinit_model=False)
			uv_fitval.append(fit.result['Parameters'][2])
			uv_fiterr.append(fit.result['Uncertainties'][2])
		else:
			print 'Please specify whether you wish to perform a Monte Carlo fit o uv, (T) or not(F)'
	    if cutout== 'T':
		immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',box=cut_reg,outfile=outputPath+label+intervalString+'_temp.image')
		immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',region='annulus['+cen_annulus+','+cen_radius+']',outfile=outputPath+label+intervalString+'_rms.image')
		#immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',box=rmsbox2,outfile=outputPath+label+intervalString+'_rms2.image')
		#immath(imagename=outputPath+label+intervalString+'.image',mode='evalexpr',expr='IM0',box=rmsbox3,outfile=outputPath+label+intervalString+'_rms3.image')
		comm_and1='rm -rf '+outputPath+label+intervalString+'.*'
		comm_and2='mv '+outputPath+label+intervalString+'_temp.image '+outputPath+label+intervalString+'.image'
		os.system(comm_and1)
		os.system(comm_and2)


uv_fitval_arr=np.array(uv_fitval)
uv_fiterr_arr=np.array(uv_fiterr)

############################################################################################
# Read values from all fitting files into lists for plotting and creating output data file
############################################################################################
# Initialise lists.
fluxDensity = []
fluxDensity2 =[]
fluxDensity3 =[]
timerange = []
fluxError = []
fluxError2 =[]
fluxError3 =[]
suffix = []
suffix2 = []
suffix3 = []

# Create flux density, flux density error, suffix and timerange vectors. In each case get values from imfit text files.
# The exact layout of the imfit outputs are somewhat unpredictable. Therefore find() is used here to reliably locate
# the relavent information.

for interval, time,interval_uv,time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
    # For some intervals the imfit may have failed, in which case the imfit text file will not contain the required information.
    # Such files will end with the string '*** FIT FAILED ***', and an if statement is used to pass over these files.
    flux_uv=[]
    fluxerr_uv=[]
    if do_monte =='T':
	FD_list=[]
	FDE_list=[]
	FD2_list=[]
	FDE2_list=[]

	fluxerr_uv=[]
	for i in range(0,len(samp_px)):
		intervalString=interval.replace(':', '.').replace('/','_')

    		imfitFile = open(outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text', 'r')
    		imfitText = imfitFile.read()
    		imfitFile.close()
    		end = imfitText[-19:-1]
		if end != '*** FIT FAILED ***':
        		# Find the position of the peak flux value.
    			if integ_fit == 'T':
    				startWord = '--- Integrated:'
        			endWord = '--- Peak:'
    				startPos = imfitText.find(startWord)
        			endPos = imfitText.find(endWord)
        			fluxString = imfitText[startPos:endPos]
    				plusMinusPos = fluxString.find('+/-')
        			unitsPos = fluxString.find('Jy')
    				plot_label_unit=')'
    				fluxD=float(fluxString[16:plusMinusPos])
				fluxE=float(fluxString[plusMinusPos+3:unitsPos-1])
				suff = str(fluxString[unitsPos-1])
				if suff == 'u':
        				fluxD=fluxD*[1.0e-6]
					fluxE=fluxE*[1.0e-6]
    				elif suff == 'm':
        				fluxD=fluxD*[1.0e-3]
					fluxE=fluxE*[1.0e-3]
    				elif suff == ' ':
        				fluxD=fluxD*[1.0]
					fluxE=fluxE*[1.0]
    				elif suff == 'n':
        				fluxD=fluxD*[1.0e-9]
					fluxE=fluxE*[1.0e-9]
				FD_list.append(fluxD)
				FDE_list.append(fluxE)
    			elif integ_fit == 'F':
        			startWord = '--- Peak:'
        			endWord = '--- Polarization:'
        			startPos = imfitText.find(startWord)
        			endPos = imfitText.find(endWord)
        			fluxString = imfitText[startPos:endPos]
    				plot_label_unit='/beam)'
        		# fluxString now contains something like "--- Peak:     432.1 +/- 4.6 uJy/beam"
       	 		# There may or may not be decimals in the values, and this tends to shift the exact
        		# positions of the values between imfits. As a workaround the '+/-' and 'Jy/b' are
        		# located as reference points.
        			plusMinusPos = fluxString.find('+/-')
        			unitsPos = fluxString.find('Jy/b')
        			fluxD=float(fluxString[10:plusMinusPos])
				fluxE=float(fluxString[plusMinusPos+3:unitsPos-1])
				suff = str(fluxString[unitsPos-1])
				if suff == 'u':
        				fluxD=fluxD*1.0e-6
					fluxE=fluxE*1.0e-6
    				elif suff == 'm':
        				fluxD=fluxD*1.0e-3
					fluxE=fluxE*1.0e-3
    				elif suff == ' ':
        				fluxD=fluxD*1.0
					fluxE=fluxE*1.0
    				elif suff == 'n':
        				fluxD=fluxD*1.0e-9
					fluxE=fluxE*1.0e-9
				FD_list.append(fluxD)
				FDE_list.append(fluxE)
			elif integ_fit == 'B':
				startWord = '--- Integrated:'
        			endWord = '--- Peak:'
    				startPos = imfitText.find(startWord)
        			endPos = imfitText.find(endWord)
        			fluxString = imfitText[startPos:endPos]
    				plusMinusPos = fluxString.find('+/-')
        			unitsPos = fluxString.find('Jy')
    				plot_label_unit=')'
    				fluxD1=float(fluxString[16:plusMinusPos])
				fluxE1=float(fluxString[plusMinusPos+3:unitsPos-1])
				suff = str(fluxString[unitsPos-1])
				if suff == 'u':
        				fluxD1=fluxD1*1.0e-6
					fluxE1=fluxE1*1.0e-6
    				elif suff == 'm':
        				fluxD1=fluxD1*1.0e-3
					fluxE1=fluxE1*1.0e-3
    				elif suff == ' ':
        				fluxD1=fluxD1*1.0
					fluxE1=fluxE1*1.0
    				elif suff == 'n':
        				fluxD1=fluxD1*1.0e-9
					fluxE1=fluxE1*1.0e-9
				startWord2 = '--- Peak:'
        			endWord2 = '--- Polarization:'
        			startPos2 = imfitText.find(startWord2)
        			endPos2 = imfitText.find(endWord2)
        			fluxString2 = imfitText[startPos2:endPos2]
    				plot_label_unit2='/beam)'
				plusMinusPos2 = fluxString2.find('+/-')
        			unitsPos2 = fluxString2.find('Jy/b')
        			fluxD2=float(fluxString2[10:plusMinusPos2])
				fluxE2=float(fluxString2[plusMinusPos2+3:unitsPos2-1])
				suff2 = str(fluxString2[unitsPos2-1])
				if suff2 == 'u':
        				fluxD2=fluxD2*1.0e-6
					fluxE2=fluxE2*1.0e-6
    				elif suff2 == 'm':
        				fluxD2=fluxD2*1.0e-3
					fluxE2=fluxE2*1.0e-3
    				elif suff2 == ' ':
        				fluxD2=fluxD2*1.0
					fluxE2=fluxE2*1.0
    				elif suff2 == 'n':
        				fluxD2=fluxD2*1.0e-9
					fluxE2=fluxE2*1.0e-9
				FD_list.append(fluxD1)
				FDE_list.append(fluxE1)
				FD2_list.append(fluxD2)
				FDE2_list.append(fluxE2)
        		#fluxE=float(fluxString[plusMinusPos+3:unitsPos-1])
        		#suffix = suffix + [str(fluxString[unitsPos-1])]
        		#timerange = timerange + [interval]
			#FD_list.append(fluxD)
			#FDE_list.append(fluxE)

	if len(FD_list) < nsim/2.:
		print '\nImage Fit failed on interval ' + interval
        	# The corresponding time intervals need to be removed from the
        	timeIntervals.remove(interval)
        	mjdTimes.remove(time)

	GaussPer=norm.cdf(1)
	median=np.percentile(FD_list,50)
	pct15=np.percentile(FD_list,(1-GaussPer)*100.0)
	pct85=np.percentile(FD_list,GaussPer*100.0)
	upp_err=pct85-median
	low_err=median-pct15
	era=(upp_err+low_err)/2.
	fluxDensity = fluxDensity + [median*1e3]
	fluxError = fluxError + [era*1e3]

	if integ_fit == 'B':
		median2=np.percentile(FD2_list,50)
		pct152=np.percentile(FD2_list,(1-GaussPer)*100.0)
		pct852=np.percentile(FD2_list,GaussPer*100.0)
		upp_err2=pct852-median2
		low_err2=median2-pct152
		era2=(upp_err2+low_err2)/2.
		fluxDensity2 = fluxDensity2 + [median2*1e3]
		fluxError2 = fluxError2 + [era2*1e3]
		suffix2 = suffix2 + ['m']

	#
	#
	timerange = timerange + [interval]
	suffix = suffix + ['m']
    		#else:
        		#print '\nFit failed on interval ' + interval
        		# The corresponding time intervals need to be removed from the
        		#timeIntervals.remove(interval)
        		#mjdTimes.remove(time)



    elif do_monte == 'F':
    	intervalString=interval.replace(':', '.').replace('/','_')

    	imfitFile = open(outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text', 'r')
    	imfitText = imfitFile.read()
    	imfitFile.close()
    	end = imfitText[-19:-1]
    	if end != '*** FIT FAILED ***':
        	# Find the position of the peak flux value.
    		if integ_fit == 'T':
    			startWord = '--- Integrated:'
        		endWord = '--- Peak:'
    			startPos = imfitText.find(startWord)
        		endPos = imfitText.find(endWord)
        		fluxString = imfitText[startPos:endPos]
    			plusMinusPos = fluxString.find('+/-')
        		unitsPos = fluxString.find('Jy')
    			plot_label_unit=')'
    			fluxDensity = fluxDensity + [float(fluxString[16:plusMinusPos])]
			fluxError = fluxError + [float(fluxString[plusMinusPos+3:unitsPos-1])]
			suffix = suffix + [str(fluxString[unitsPos-1])]
    		elif integ_fit == 'F':
        		startWord = '--- Peak:'
        		endWord = '--- Polarization:'
        		startPos = imfitText.find(startWord)
        		endPos = imfitText.find(endWord)
        		fluxString = imfitText[startPos:endPos]
    			plot_label_unit='/beam)'
        	# fluxString now contains something like "--- Peak:     432.1 +/- 4.6 uJy/beam"
       	 	# There may or may not be decimals in the values, and this tends to shift the exact
        	# positions of the values between imfits. As a workaround the '+/-' and 'Jy/b' are
        	# located as reference points.
        		plusMinusPos = fluxString.find('+/-')
        		unitsPos = fluxString.find('Jy/b')
        		fluxDensity = fluxDensity + [float(fluxString[10:plusMinusPos])]
			fluxError = fluxError + [float(fluxString[plusMinusPos+3:unitsPos-1])]
			suffix = suffix + [str(fluxString[unitsPos-1])]
		elif integ_fit == 'B':
				startWord = '--- Integrated:'
        			endWord = '--- Peak:'
    				startPos = imfitText.find(startWord)
        			endPos = imfitText.find(endWord)
        			fluxString = imfitText[startPos:endPos]
    				plusMinusPos = fluxString.find('+/-')
        			unitsPos = fluxString.find('Jy')
    				plot_label_unit=')'
    				fluxDensity= fluxDensity +[float(fluxString[16:plusMinusPos])]
				fluxError=fluxError +[float(fluxString[plusMinusPos+3:unitsPos-1])]
				suffix = suffix + [str(fluxString[unitsPos-1])]
				startWord2 = '--- Peak:'
        			endWord2 = '--- Polarization:'
        			startPos2 = imfitText.find(startWord2)
        			endPos2 = imfitText.find(endWord2)
        			fluxString2 = imfitText[startPos2:endPos2]
    				plot_label_unit2='/beam)'
				plusMinusPos2 = fluxString2.find('+/-')
        			unitsPos2 = fluxString2.find('Jy/b')
        			fluxDensity2=fluxDensity2+[float(fluxString2[10:plusMinusPos2])]
				fluxError2=fluxError2+[float(fluxString2[plusMinusPos2+3:unitsPos2-1])]
				suffix2 = suffix2 + [str(fluxString2[unitsPos2-1])]
        	#fluxError = fluxError + [float(fluxString[plusMinusPos+3:unitsPos-1])]
        	#suffix = suffix + [str(fluxString[unitsPos-1])]
        	timerange = timerange + [interval]
    	else:
        	print '\nFit failed on interval ' + interval
        	# The corresponding time intervals need to be removed from the
		fluxError_real.remove(fluxError_real[timeIntervals.index(interval)])
        	timeIntervals.remove(interval)
        	mjdTimes.remove(time)

    else:
	print 'Please specify whether you wish to perform a Monte Carlo fit (T) or not(F)'
    GaussPer=norm.cdf(1)
    if uv_fit == 'T'and do_monte_uv=='T':
		median3=np.percentile(uv_fitval_arr[np.where(np.array(timeIntervals)==interval)[0]][0],50)
		pct153=np.percentile(uv_fitval_arr[np.where(np.array(timeIntervals)==interval)[0]][0],(1-GaussPer)*100.0)
		pct853=np.percentile(uv_fitval_arr[np.where(np.array(timeIntervals)==interval)[0]][0],GaussPer*100.0)
		upp_err3=pct853-median3
		low_err3=median3-pct153
		era3=(upp_err3+low_err3)/2.
		fluxDensity3 = fluxDensity3 + [median3]
		fluxError3 = fluxError3 + [era3]
		suffix3 = suffix3 + [' ']

		if len(uv_fitval[np.where(np.array(timeIntervals)==interval)[0]]) < nsim/2.:
			print '\nUV Fit failed on interval ' + interval
        		# The corresponding time intervals need to be removed from the
        		#timeIntervals_uv.remove(interval_uv)
        		#mjdTimes_uv.remove(time_uv)

    if uv_fit == 'T'and do_monte_uv=='F':
		flux_uv1=uv_fitval_arr[np.where(np.array(timeIntervals)==interval)[0]][0]
		fluxerr_uv1=uv_fiterr_arr[np.where(np.array(timeIntervals)==interval)[0]][0]
		fluxDensity3 = fluxDensity3 + [flux_uv1]
		fluxError3 = fluxError3 + [fluxerr_uv1]
		suffix3 = suffix3 + [' ']

########################################################################
#Print results to stdout for user--> probably should get rid of this,
#maybe just a 'script successful' message?
#########################################################################
#print image plane fitting results for user
if integ_fit == 'B':
	print '\nFlux desities extracted from imfit (integ):'
	print fluxDensity
	print '\nFlux errors extracted from imfit (integ):'
	print fluxError
	print '\nUnits extracted from imfit (integ):'
	print suffix
	print '\nFlux desities extracted from imfit (peak):'
	print fluxDensity2
	print '\nFlux errors extracted from imfit (peak):'
	print fluxError2
	print '\nUnits extracted from imfit (peak):'
	print suffix2
	print '\nRealistic Error from regions (Jy):'
	print fluxError_real
else:
	print '\nFlux desities extracted from imfit:'
	print fluxDensity
	print '\nFlux errors extracted from imfit:'
	print fluxError
	print '\nUnits extracted from imfit:'
	print suffix
	print '\nRealistic Error from regions (Jy):'
	print fluxError_real

#print uv fitting results for user, if requested
if uv_fit == 'T':
	print '\nFlux desities extracted from uvmodelfit:'
	print fluxDensity3
	print '\nFlux errors extracted from uvmodelfit:'
	print fluxError3
	print '\nUnits extracted from uvmodelfit:'
	print suffix3

########################################################################
#Make sure verything is in the right units
#########################################################################
# The flux values will be plotted in uJy/bm, however the imfit outputs will not necessarily be in these
# units. To get around this the units for each value are recorded in the suffix list created in the above
# for loop. The following loop will use this list to create a multiplier list, which will be used to
# convert all values to Jy/bm. They can then all be scaled to uJy/bm ready for plotting.
multiplier = []
multiplier2 = []
multiplier3 = []
for x in suffix:
    if x == 'u':
        multiplier = multiplier + [1.0e-6]
    elif x == 'm':
        multiplier = multiplier + [1.0e-3]
    elif x == ' ':
        multiplier = multiplier + [1.0]
    elif x == 'n':
        multiplier = multiplier + [1.0e-9]
print multiplier
if integ_fit=='B':
	for x in suffix2:
    		if x == 'u':
        		multiplier2 = multiplier2 + [1.0e-6]
    		elif x == 'm':
        		multiplier2 = multiplier2 + [1.0e-3]
    		elif x == ' ':
        		multiplier2 = multiplier2 + [1.0]
    		elif x == 'n':
			multiplier2 = multiplier2 + [1.0e-9]
	print 'peak', multiplier2
if uv_fit=='T':
	for x in suffix3:
    		if x == 'u':
        		multiplier3 = multiplier3 + [1.0e-6]
    		elif x == 'm':
        		multiplier3 = multiplier3 + [1.0e-3]
    		elif x == ' ':
        		multiplier3 = multiplier3 + [1.0]
    		elif x == 'n':
			multiplier3 = multiplier3 + [1.0e-9]
	print 'uv peak', multiplier3

# Multiplier() should now contain a list of multiplication factors ready to be applied to the flux data. The
#end result is in Jy
length = range(len(fluxDensity))
for n in length:
    if integ_fit == 'B':
	fluxDensity2[n] = fluxDensity2[n] * multiplier2[n]
    	fluxError2[n] = fluxError2[n] * multiplier2[n]
    if uv_fit == 'T':
	fluxDensity3[n] = fluxDensity3[n] * multiplier3[n]
    	fluxError3[n] = fluxError3[n] * multiplier3[n]
    fluxDensity[n] = fluxDensity[n] * multiplier[n]
    fluxError[n] = fluxError[n] * multiplier[n]

# The flux values will be plotted in a user specified unit, so they are scaled here,
#along with error values.
for k in length:
    if integ_fit == 'B':
	fluxDensity2[k] = fluxDensity2[k]*lc_scale_factor
    	fluxError2[k] = fluxError2[k]*lc_scale_factor
    if uv_fit == 'T':
	fluxDensity3[k] = fluxDensity3[k]*lc_scale_factor
    	fluxError3[k] = fluxError3[k]*lc_scale_factor
    fluxDensity[k] = fluxDensity[k]*lc_scale_factor
    fluxError[k] = fluxError[k]*lc_scale_factor
    fluxError_real[k]=fluxError_real[k]*lc_scale_factor


########################################################################
#Write results to data file
#########################################################################
data = open(dataPath, 'w')
for i in range(0,len(fluxDensity)):
	if integ_fit == 'B':
		data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|peak flux|peak imfit err|rms error'))
		data.write('{0} {1} {2} {3} {4} {5}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity2[i],fluxError2[i],fluxError_real))
		if uv_fit == 'T':
			data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|peak flux|peak imfit err|uv flux|uverr|rms error'))
			data.write('{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity2[i],fluxError2[i],fluxDensity3[i],fluxError3[i],fluxError_real))
	elif integ_fit =='T':
		data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|rms error'))
		data.write('{0} {1} {2} {3}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxError_real))
		if uv_fit=='T':
			data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|uv flux|uverr|rms error'))
			data.write('{0} {1} {2} {3} {4} {5}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity3[i],fluxError3[i],fluxError_real))
	elif integ_fit =='F':
		data.write('{0}\n'.format('#MJD|peak flux|peak imfit err|rms error'))
		data.write('{0} {1} {2} {3}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxError_real))
		if uv_fit=='T':
			data.write('{0}\n'.format('#MJD|peak flux|peak imfit err|uv flux|uverr|rms error'))
			data.write('{0} {1} {2} {3} {4} {5}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity3[i],fluxError3[i],fluxError_real))
data.close()

########################################################################
#Plot Lightcurves
#########################################################################
minutesElapsed=[]
secondsElapsed=[]
hoursElapsed=[]
for i in range(len(mjdTimes)):
    hoursElapsed.append((mjdTimes[i]-mjdTimes[0])*24+intervalSizeS/(60.0*60.0*2.0))
    minutesElapsed.append((mjdTimes[i]-mjdTimes[0])*24*60+intervalSizeS/(60.0*2.0))
    secondsElapsed.append((mjdTimes[i]-mjdTimes[0])*24*60*60+intervalSizeS/2.0)

if lc_scale_time=='M':
	Elapsed=minutesElapsed
elif lc_scale_time=='S':
	Elapsed=secondsElapsed
elif lc_scale_time=='H':
	Elapsed=hoursElapsed

if integ_fit == 'B':
	fig1=pp.figure()
	pp.errorbar(Elapsed, fluxDensity, yerr=fluxError_real, fmt='ro',)
	pp.xlabel('Time since start of observation (mins)')
	y_label_name='Flux Density (mJy'+plot_label_unit+')'
	pp.ylabel(y_label_name)
	pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	savestring = path_dir+'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_check_lc_integ.eps'
	pp.savefig(savestring)
	print savestring, ' is saved'
	fig2=pp.figure()
	pp.errorbar(Elapsed, fluxDensity2, yerr=fluxError_real, fmt='ro',)
	pp.xlabel('Time since start of observation (mins)')
	y_label_name2='Flux Density (mJy'+plot_label_unit2+')'
	pp.ylabel(y_label_name2)
	pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	savestring2 = path_dir+'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_check_lc_peak.eps'
	pp.savefig(savestring2)
	print savestring2, ' is saved'
	if uv_fit=='T':
		fig3=pp.figure()
		pp.errorbar(Elapsed, fluxDensity3, yerr=fluxError3, fmt='ro',)
		pp.xlabel('Time since start of observation (mins)')
		y_label_name='Flux Density (mJy/beam)'
		pp.ylabel(y_label_name)
		pp.title('Flux Density vs Time. '+target+' '+refFrequency)
		pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
		savestring = path_dir+'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_check_lc_uv.eps'
		pp.savefig(savestring)
		print savestring, ' is saved'

else:
	pp.errorbar(Elapsed, fluxDensity, yerr=fluxError_real, fmt='ro',)
	pp.xlabel('Time since start of observation (mins)')
	y_label_name='Flux Density (mJy'+plot_label_unit+')'
	pp.ylabel(y_label_name)
	pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	savestring = path_dir+'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_check_lc.eps'
	pp.savefig(savestring)
	print savestring, ' is saved'
	if uv_fit=='T':
		fig3=pp.figure()
		pp.errorbar(Elapsed, fluxDensity3, yerr=fluxError3, fmt='ro',)
		pp.xlabel('Time since start of observation (mins)')
		y_label_name='Flux Density (mJy/beam)'
		pp.ylabel(y_label_name)
		pp.title('Flux Density vs Time. '+target+' '+refFrequency)
		pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
		savestring = path_dir+'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_check_lc_uv.eps'
		pp.savefig(savestring)
		print savestring, ' is saved'
#remove temp files and .last/.log files created by CASA/ipython
os.system('rm -rf *.last')
os.system('rm -rf *.log')
os.system('rm -rf tempfile.txt tempfile2.txt')
