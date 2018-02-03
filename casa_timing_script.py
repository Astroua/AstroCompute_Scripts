##################################
#CASA Timing Script
##################################
''' CASA script to produce high time res light curves from calibrated continuum MS
INPUT: Calibrated and Split MS, parameter file (param_file below)
OUTPUT: (1) Lightcurve plot--> pdf saved in [path_dir]/data_products
        (2) Lightcurve data file--> saved in [path_dir]/data_products
        (3) (optional) Lomb-Scarge periodogram plot and basic variability properties file --> saved in [path_dir]/data_products
NOTES: - This script is theoretically compatible with any data that can be imported into a CASA MS,
         but has only been tested with VLA, SMA, and NOEMA data.

Written by: C. Gough (original version), additions and updates by A. Tetarenko & E. Koch
Last Updated: Feb 2018

TO RUN SCRIPT-->casa -c casa_timing_script.py [path_to_param_file] [path_dir] [path_to_repo]
Uncomment at line 557 if you dont want time-bins printed to screen
Uncomment at line 1039 if you dont want results printed to screen.

NOTE: path_dir is path to input/output directory of your choice
-MS's need to be in path_dir/data,
-all output goes to path_dir/data_products (this directory is created by script)
-Make sure to including trailing / in [path_to_repo_dir] & [path_dir]!!!
-Remember to put Aegean directory location in python path if using object detection.
'''

#MODULES USED:
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
import astroML.time_series
import imp
import warnings
warnings.filterwarnings('ignore')

import sys
path_ac = sys.argv[-1]
if os.path.isdir(path_ac):
    sys.path.append(path_ac)
else:
    raise Exception('Please set the path to the AstroCompute_Scripts directory correctly.')

from utils import (convert_param_format,initial_clean,run_aegean,
                   var_analysis,lomb_scargle,chi2_calc,errf)
from Aegean_ObjDet import objdet

##################################
#Setup
##################################
path_dir = sys.argv[-2]
param_file = sys.argv[-3]

if not os.path.isdir(os.path.join(path_dir, 'data_products/')):
    os.mkdir(os.path.join(path_dir, "data_products"))
if not os.path.isdir(os.path.join(path_dir, 'data/')):
    os.mkdir(os.path.join(path_dir, "data"))
    raise Exception(path_dir+'data/ created. Please move your MS into that directory.')

#get input parameters from file
#from utils import load_json
#data_params = load_json(param_file)
##to convert txt param file directly to dictionary do this instead:
data_params = convert_param_format(param_file, to="dict")

##################################
#Reading in Parameters
##################################
print 'Reading in parameters...\n'
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
intervalSizeH = int(data_params["intervalSizeH"])
intervalSizeM = int(data_params["intervalSizeM"])
intervalSizeS = float(data_params["intervalSizeS"])
frac,whole=m.modf(float(intervalSizeS))
intervalSizemicro = int(frac*(1e6))
intervalSizeSec = int(whole)
# Name of visibliity
visibility = os.path.join(path_dir, 'data/' + data_params["visibility"])
# do you want to do uv fitting (T or F) as well?
uv_fit = data_params["uv_fit"]
#make copy for uv fitting
visibility_uv= visibility.rstrip('.ms') + '_uv.ms'
if uv_fit=='T':
    if not os.path.isdir(visibility_uv):
    	os.system('cp -r '+visibility+' '+visibility_uv)

''' VARIABILITY ANALYSIS'''
#Do you want a basic variability analysis?
var_anal=data_params["var_anal"]
power_spec=data_params["power_spec"]
#Variability file name
dataPathVar = \
    os.path.join(path_dir,'data_products/varfile_'+target+ '_' + obsDate +'_'+refFrequency +
                 '_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min' +
                 str(intervalSizeS)+'sec.txt')
#periodogram name
labelP = \
    os.path.join(path_dir,'data_products/periodogram_'+target+ '_' + obsDate +'_'+refFrequency +
                 '_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min' +
                 str(intervalSizeS)+'sec.pdf')

'''DIRECTORY AND FILE NAME PARAMETERS'''
# Set path to directory where all output from this script is saved.
outputPath = os.path.join(path_dir, 'data_products/images_'+target+ '_' + obsDate +'_'+refFrequency+'_'+\
str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec/')

# dataPath contains the path and filename in which data file will be saved.
# This script can be run on several epochs of data from the same observation without changing this path.
# In this case the data file will be appended each time.
dataPath = \
    os.path.join(path_dir,'data_products/datafile_'+target+ '_' + obsDate +'_'+refFrequency +
                 '_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min' +
                 str(intervalSizeS)+'sec.txt')

#make output directory (within data_products directory)--> done in initial clean script too
if not os.path.isdir(os.path.join(path_dir,
                                  'data_products/images_'+target+ '_' + obsDate +'_' +
                                  refFrequency+'_'+str(intervalSizeH) +
                                  'hours'+str(intervalSizeM)+'min' +
                                  str(intervalSizeS)+'sec')):
	mkdir_string='mkdir '+ os.path.join(path_dir, 'data_products/images_'+target+ '_' + obsDate +'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec')
	mkdir_perm1='sudo chown ubuntu '+ os.path.join(path_dir, 'data_products/images_'+target+ '_' + obsDate +'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec')
	mkdir_perm2='sudo chmod -R u+rwx '+ os.path.join(path_dir, 'data_products/images_'+target+ '_' + obsDate +'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec')
	os.system(mkdir_string)
	# os.system(mkdir_perm1)
	# os.system(mkdir_perm2)

'''FLAGS'''
# flag to run object detection
runObj = data_params["runObj"]
# Is the data set large enough that you only want to save a cutout? Always yes now!
cutout = 'T'
# Clean can be toggled on/off here (T/F/U).
runClean = data_params["runClean"]
# optimize cleaning
opt_clean = 'F'#data_params["opt_clean"]
# do you want to fix parameters in fits from full data set fit? (T of F)
fix_pos = data_params["fix_pos"]
# if fixed parameters do you want to mc sample the fixed parameters (T) or take the best fit (F)?
do_monte = data_params["do_monte"]
# do you want peak (mJy/beam; F) or integrated (mJy; T) flux, or both(B) in lightcurve file?
integ_fit = data_params["integ_fit"]
# fix position in UV
uv_fix = data_params["uv_fix"]
# If runClean=F set fit_cutout='T' to make sure all coordinates correct
if runClean=='F':
	fit_cutout = 'T'
else:
	fit_cutout = 'F'
# define start and end time yourself
def_times = data_params["def_times"]

'''CLEAN PARAMETERS'''
# The clean command (line 530) should be inspected closely to ensure all arguments are appropriate before
# running script on a new data set.
# The following arguments will be passed to casa's clean, imfit or imstat functions:
imageSize = [int(data_params["imageSize"])] * 2
if opt_clean == 'T':
	if not is_power2(imageSize[0]):
		print 'Clean will run faster if image size is 2^n'
		imageSize = \
            [int(pow(2, m.ceil(np.log(imageSize[0])/np.log(2)))),
             int(pow(2, m.ceil(np.log(imageSize[0])/np.log(2))))]
        	print 'imagesize is now set to ', imageSize
    	print 'imagesize remains at ', imageSize, 'due to user request'
numberIters = int(data_params["numberIters"])
cellSize = [data_params["cellSize"]] * 2
taylorTerms = int(data_params["taylorTerms"])
myStokes = data_params["myStokes"]
thre = data_params["thre"]
robust = float(data_params["robust"])
weighting= data_params["weighting"]
#put threshold in Jy for image convergence test below
thre_unit=re.findall("[a-zA-Z]+", thre)
if thre_unit == 'uJy':
    thre_num=float(re.findall(r"[-+]?\d*\.\d+|\d+",thre)[0])/1e6
elif thre_unit == 'mJy':
    thre_num=float(re.findall(r"[-+]?\d*\.\d+|\d+",thre)[0])/1e3
elif thre_unit == 'nJy':
    thre_num=float(re.findall(r"[-+]?\d*\.\d+|\d+",thre)[0])/1e3
else:
    thre_num=float(re.findall(r"[-+]?\d*\.\d+|\d+",thre)[0])/1.
spw_choice = data_params["spw_choice"]
# If an outlier file is to be used in the clean, set outlierFile to the filename (path inluded). thre will
# also need to be set if outlier file is to be used.
# If not, set outlierFile to ''.
outlierFile = data_params["outlierFile"]


'''OBJECT DETECTION AND SELECTION PARAMETERS'''
mask_option=data_params["mask_option"]
if runObj == 'T':
	tables = outputPath+label+'whole_dataset_objdet_comp.tab'
	out_file0=outputPath+label+'whole_dataset_aegean.txt'
	catalog_input_name=outputPath+label+'whole_dataset_objdet.tab'
	tele=data_params["tele"]
	lat=data_params["lat"]
	seed=int(data_params["seed"])
	flood=int(data_params["flood"])
	fits_file=outputPath+label+'whole_dataset.fits'
	initial_clean(visibility,outputPath,label,imageSize,cellSize,spw_choice,taylorTerms,numberIters,thre,robust,weighting)
	src_l,ra_l,dec_l,maj_l,min_l,pos_l=objdet(tele,lat,out_file0,fits_file,seed,flood,tables,catalog_input_name,cellSize[0])
	if len(src_l)==0:
		raise Exception('No sources detected.')
	elif len(src_l)==1:
		ind=1
	else:
		ind=raw_input('Please enter target source ID--> ')
	# target position-->Take bounding ellipse from Aegean and convert to minimum bounding box in pixels for
	#use with rest of script
	mask_file_OD=open(outputPath+'aegean_mask.txt','w')
	mask_file_OD.write('{0}\n'.format('#CRTFv0'))
	for i in range(0,len(src_l)):
		pos=au.findRADec(outputPath+label+'whole_dataset.image',ra_l[int(i)-1]+' '+dec_l[int(i)-1])
		bbox_halfwidth=np.sqrt((min_l[int(i)-1]*np.cos(pos_l[int(i)-1]))**2+(min_l[int(i)-1]*np.sin(pos_l[int(i)-1]))**2)+3
		bbox_halfheight=np.sqrt((maj_l[int(i)-1]*np.cos(pos_l[int(i)-1]+(np.pi/2.)))**2+(maj_l[int(i)-1]*np.sin(pos_l[int(i)-1]+(np.pi/2.)))**2)+3
		Box = str(pos[0]-bbox_halfwidth)+','+ str(pos[1]-bbox_halfheight)+','+str(pos[0]+bbox_halfwidth)+','+ str(pos[1]+bbox_halfheight)
		mask_file_OD.write('{0}\n'.format('box[['+Box.split(',')[0]+'pix,'+Box.split(',')[1]+'pix], ['+Box.split(',')[2]+'pix,'+Box.split(',')[3]+'pix]]'))
	mask_file_OD.close()
	tar_pos=au.findRADec(outputPath+label+'whole_dataset.image',ra_l[int(ind)-1]+' '+dec_l[int(ind)-1])
	bbox_halfwidth=np.sqrt((min_l[int(ind)-1]*np.cos(pos_l[int(ind)-1]))**2+(min_l[int(ind)-1]*np.sin(pos_l[int(ind)-1]))**2)+3
	bbox_halfheight=np.sqrt((maj_l[int(ind)-1]*np.cos(pos_l[int(ind)-1]+(np.pi/2.)))**2+(maj_l[int(ind)-1]*np.sin(pos_l[int(ind)-1]+(np.pi/2.)))**2)+3
	targetBox = str(tar_pos[0]-bbox_halfwidth)+','+ str(tar_pos[1]-bbox_halfheight)+','+str(tar_pos[0]+bbox_halfwidth)+','+ str(tar_pos[1]+bbox_halfheight)
	mask_option='aegean'
elif runObj == 'F':
    # input target box in pixels if not running object detection
    targetBox = data_params["targetBox"] # #'2982,2937,2997,2947'
else:
    raise ValueError("runObj must be 'T' or 'F'. Value given is ", runObj)

#mask for clean based on target box or mask file for complicated fields
if mask_option == 'box':
	maskPath = 'box [['+targetBox.split(',')[0]+'pix,'+targetBox.split(',')[1]+'pix],['+targetBox.split(',')[2]+\
    'pix,'+targetBox.split(',')[3]+'pix]]'
elif mask_option == 'file':
	maskPath = os.path.join(path_dir, 'data/'+data_params["mask_file"])
elif mask_option == 'aegean':
	maskPath=os.path.join(outputPath, 'aegean_mask.txt')
else:
	raise ValueError("mask_option must be 'box' or 'file'. Value given is ", mask_option)

'''IMAGE PRODUCT PARAMETERS: CUTOUT AND RMS/ERROR IN IMAGE PLANE HANDLING'''
#define rms boxes for realistic error calculation:
#pix_shift_cutout is how many pixels past target box you want in cutout images rms calculation
#annulus_rad_inner/outer is inner and outer radius in pixels from source
pix_shift_cutout=30
annulus_rad_inner=20
annulus_rad_outer=30
cut_reg=str(float(targetBox.split(',')[0])-float(pix_shift_cutout))+','+str(float(targetBox.split(',')[1])-float(pix_shift_cutout))+\
','+str(float(targetBox.split(',')[2])+float(pix_shift_cutout))+','+str(float(targetBox.split(',')[3])+float(pix_shift_cutout))#'2962,2917,3017,2967'
x_sizel=float(targetBox.split(',')[0])
x_sizeu=float(targetBox.split(',')[2])
y_sizel=float(targetBox.split(',')[1])
y_sizeu=float(targetBox.split(',')[3])
#annulus region parameters
cen_annulus='['+str(((x_sizeu-x_sizel)/2.)+x_sizel)+'pix,'+str(((y_sizeu-y_sizel)/2.)+y_sizel)+'pix]'
cen_radius='['+str(annulus_rad_inner)+'pix,'+str(annulus_rad_outer)+'pix]'
#what unit do you want light curve
lc_scale_unit=data_params["lc_scale_unit"]
if lc_scale_unit=='m':
	lc_scale_factor=1.0e03
elif lc_scale_unit=='u':
	lc_scale_factor=1.0e06
elif lc_scale_unit=='n':
    lc_scale_factor=1.0e09
elif lc_scale_unit=='':
	lc_scale_factor=1.0
else:
    lc_scale_factor=1.0
#what time scale do you want in light curve
lc_scale_time=data_params["lc_scale_time"]
#if remainder of obstime/interval > this # then it is included in time intervals (in minutes)
rem_int=5.



'''IMAGE FITTING PARAMETERS'''
if do_monte == 'T':
#add appropriate label to final lightcurve file name
	lab='_mc_'
else:
	lab='_bestfit_'
# number of MC simulations if MC position sampling chosen
nsim=100
#if fixing parameters, what parameters do you want to fix in fits? f= peak flux, x=peak x pos, y=peak y pos,
#a=major axis (arcsec), b=minor axis (arcsec), p=position angle (deg)
#a, b, p "convolved with beam" values not "deconvolved" values!!!
#format is [value,error] from fit of full data set
if fix_pos == 'T':
	par_fix='xy'#data_params["par_fix"]
	if runObj == 'F':
		print 'Cleaning Full Data Set-->'
		os.system('rm -rf '+outputPath+label+'whole_dataset.*')
		initial_clean(visibility,outputPath,label,imageSize,cellSize,spw_choice,taylorTerms,numberIters,thre,robust,weighting)
	print 'Fitting full data set in Image Plane-->'
	os.system('rm -rf '+outputPath+label+'whole_dataset.txt')
	full_fit=imfit(imagename=outputPath+label+'whole_dataset.image',box=targetBox,logfile=outputPath+label+'whole_datasetfit.txt')
	imfitFilefull=open(outputPath+label+'whole_datasetfit.txt','r')
	for line in imfitFilefull:
		if ('--- ra:' in line) & ('pixels' in line):
			ra_string=line
		if ('--- dec:' in line) & ('pixels' in line):
			dec_string=line
	pmra=ra_string.find('+/-')
	pmdec=dec_string.find('+/-')
	pixra=ra_string.find('pixels')
	pixdec=dec_string.find('pixels')
	rara=ra_string.find('ra:')
	decdec=dec_string.find('dec:')
	imfitFilefull.close()

	peak_x=[float(ra_string[rara+4:pmra]),float(ra_string[pmra+3:pixra])]#[2988.63,0.02]
	peak_y=[float(dec_string[decdec+4:pmdec-1]),float(dec_string[pmdec+3:pixdec-1])]#[2942.09,0.01]
	if fit_cutout=='T':
		peak_x=[peak_x[0]-float(cut_reg.split(',')[0]),peak_x[1]]
		peak_y=[peak_y[0]-float(cut_reg.split(',')[1]),peak_y[1]]
	b_maj=[full_fit['results']['component0']['shape']['majoraxis']['value'],\
	full_fit['results']['component0']['shape']['majoraxiserror']['value']]#[0.154,0.001]
	b_min=[full_fit['results']['component0']['shape']['minoraxis']['value'],\
	full_fit['results']['component0']['shape']['minoraxiserror']['value']]#[0.099,0.0005]
	pos_ang=[full_fit['results']['component0']['shape']['positionangle']['value'],\
	full_fit['results']['component0']['shape']['positionangleerror']['value']]#[67.41,0.42]

'''REFITTTING WITHOUT CLEANING PARAMETERS'''
#make sure rms boxes are set to local region.
if fit_cutout == 'T':
	targetBox = str(float(targetBox.split(',')[0])-float(cut_reg.split(',')[0]))+','+\
    str(float(targetBox.split(',')[1])-float(cut_reg.split(',')[1]))+','+\
    str(float(targetBox.split(',')[2])-float(cut_reg.split(',')[0]))+','+\
    str(float(targetBox.split(',')[3])-float(cut_reg.split(',')[1]))
	cx_sizel=float(targetBox.split(',')[0])
	cx_sizeu=float(targetBox.split(',')[2])
	cy_sizel=float(targetBox.split(',')[1])
	cy_sizeu=float(targetBox.split(',')[3])
	cen_annulus='['+str(((cx_sizeu-cx_sizel)/2.)+cx_sizel)+'pix,'+str(((cy_sizeu-cy_sizel)/2.)+cy_sizel)+'pix]'
	cen_radius='['+str(annulus_rad_inner)+'pix,'+str(annulus_rad_outer)+'pix]'

'''UV FITTING PARAMETERS'''
if runClean=='U':
    uv_fit='T'
#only point sources right now
comp_uv='delta'
#for VLA 'I,Q,U,V', for SMA 'LL' (from listobs)
stokes_param=myStokes#data_params["stokes_param"]
if uv_fit=='T':
	print 'Fitting full data set in UV Plane-->'
	combuv=visibility_uv.strip('.ms')
	mstransform(vis=visibility_uv,outputvis=combuv+'_mstransform.ms',combinespws=True,spw='',datacolumn='data')
	fitfulluv=uvm.uvmultifit(vis=combuv+'_mstransform.ms', spw=spw_choice, column = "data", uniform=False, model=[comp_uv],stokes = stokes_param, var=['p[0],p[1],p[2]'],outfile = outputPath+label+'whole_dataset_uv.txt', OneFitPerChannel=False ,cov_return=False,finetune=False, method="levenberg")
	if uv_fix=='F':
		uv_var='p[0],p[1],p[2]'#'2.8194e-02,8.5502e-03,p[2]'
	elif uv_fix=='T':
		uv_var=str(fitfulluv.result['Parameters'][0])+','+str(fitfulluv.result['Parameters'][1])+',p[2]'#'2.8194e-02,8.5502e-03,p[2]'
	src_uv_init=str(fitfulluv.result['Parameters'][0])+','+str(fitfulluv.result['Parameters'][1])+','+\
	str(fitfulluv.result['Parameters'][2])#[2.8194e-02,8.5502e-03 , 1.3508e-01]
	src_uv_err=str(fitfulluv.result['Uncertainties'][0])+','+str(fitfulluv.result['Uncertainties'][1])+','+\
	str(fitfulluv.result['Uncertainties'][2])#[4.7722e-05 , 3.7205e-05, 1.1192e-04]
##################################


##################################
#Calculating time bins
##################################
print 'Calculating time-bins...\n'
# Run listobs and store as a text file.
if os.path.isfile(outputPath + label + 'listobs.text'):
	print 'listobs already exits.'
else:
	listobs(vis=visibility, listfile=outputPath + label + 'listobs.text')

# Get time range of observation from listobs and generate timeIntervals vector.
# Extract start and end times and start date from listobs.
file_listobs=open(outputPath + label + 'listobs.text','r')
readit=file_listobs.read()
matched_lines = [line for line in readit.splitlines() if "Observed from" in line]
listobsLine7=matched_lines[0]
file_listobs.close()

print 'Target Source', listobsLine7

startTimeH = int(listobsLine7[31:33])
startTimeM = int(listobsLine7[34:36])
startTimeS = float(listobsLine7[37:39])
endTimeH = int(listobsLine7[61:63])
endTimeM = int(listobsLine7[64:66])
endTimeS = float(listobsLine7[67:69])
startDate = listobsLine7[19:30]
endDate = listobsLine7[49:60]

# The data will be plotted against MJD time, so we need startDate in this format also.
# Convert month to integer.
year = int(startDate[7:11])
month = startDate[3:6]
day = int(startDate[0:2])
months = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
          'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
monthInt = months[month]
startDateMJD = gcal2jd(year, monthInt, day)
yeare = int(endDate[7:11])
monthe = endDate[3:6]
daye = int(endDate[0:2])
monthInte = months[monthe]
endDateMJD = gcal2jd(yeare, monthInte, daye)

####################################################################################################
# Optional: Define start and end times manually (i.e., if only want subset of observation )
#
if def_times=='T':
	day = int(data_params["startD"])
	monthInt = int(data_params["startM"])
	year = int(data_params["startY"])
	startTimeH = int(data_params["startTimeH"])
	startTimeM = int(data_params["startTimeM"])
	startTimeS = float(data_params["startTimeS"])
	daye = int(data_params["endD"])
	monthInte = int(data_params["endM"])
	yeare = int(data_params["endY"])
	endTimeH = int(data_params["endTimeH"])
	endTimeM = int(data_params["endTimeM"])
	endTimeS = float(data_params["endTimeS"])
	startDateMJD = gcal2jd(year, monthInt, day)
	endDateMJD = gcal2jd(yeare, monthInte, daye)
####################################################################################################

#This script is capable of microsecond resolution. Input is time interval in float seconds,
#so converted to microseconds here for datetime and timedelta objects below.
frac1,whole1=m.modf(float(startTimeS))
startTimeMicro = int(frac1*(1e6))
startTimeSec = int(whole1)
frac2,whole2=m.modf(float(endTimeS))
endTimeMicro = int(frac2*(1e6))
endTimeSec = int(whole2)

# We require a list of time interval strings to pass to the clean function.
# In order to generate this list we need to perform some basic arithmetic on time values.
# The datetime module achieves this using datetime and timedelta objects.

# First the start and end times are converted to datetime objects.
startTime = datetime(year,monthInt,day,startTimeH,startTimeM,startTimeSec,startTimeMicro)
endTime = datetime(yeare,monthInte,daye,endTimeH,endTimeM,endTimeSec,endTimeMicro)
endTimeHMS = time(endTimeH,endTimeM,endTimeSec,endTimeMicro)
#
# Use timedelta to calculate duration.
observationDurationDelta = endTime-startTime
#
# Duration is printed in iso format as a string for the user.
print '\nTotal observation time is ' + str(observationDurationDelta)
#
# The interval length is converted to a timedelta object
intervalSize = time(intervalSizeH, intervalSizeM, intervalSizeSec,intervalSizemicro)
intervalSizeDelta = timedelta(hours=intervalSizeH, minutes=intervalSizeM, seconds=intervalSizeSec,\
    microseconds=intervalSizemicro)
#
# Calculate number of intervals. This is done by converting the duration and interval size into seconds
durationSeconds = (observationDurationDelta.days)*(3600.*24.) + (observationDurationDelta.seconds)
intervalSeconds = (intervalSize.hour)*3600 + (intervalSize.minute)*60 + (intervalSize.second)+\
(intervalSize.microsecond*1e-6)
if intervalSeconds > durationSeconds:
    raise Exception("The observation duration (", str(durationSeconds),
                    "sec) is less than the given interval (",
                    str(intervalSeconds), "sec). Decrease the interval time.")
numIntervals = int(durationSeconds / intervalSeconds)
print 'The number of time bins is: ',numIntervals

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
    date_conv=jd2gcal(startDateMJD[0],startDateMJD[1])
    month_0={1: '01', 2: '02', 3: '03', 4: '04', 5: '05', 6: '06', 7: '07', 8: '08',\
         9: '09', 10: '10', 11: '11', 12: '12'}[date_conv[1]]
    date_conv2=jd2gcal((startDateMJD[0]),(startDateMJD[1]+endTimeDelta.days))
    month_2={1: '01', 2: '02', 3: '03', 4: '04', 5: '05', 6: '06', 7: '07',\
         8: '08', 9: '09', 10: '10', 11: '11', 12: '12'}[date_conv2[1]]
    if float(endTimeDelta.days) == 0:
    	mjdTimes[element] = startDateMJD[1] + long((datetime.min+endTimeDelta).time().hour)/24.0 + \
        long((datetime.min+endTimeDelta).time().minute)/60.0/24.0 + \
        long((datetime.min+endTimeDelta).time().second)/60.0/60.0/24.0+ \
        long((datetime.min+endTimeDelta).time().microsecond)/60.0/60.0/24.0/(1.0e6)
    else:
    	mjdTimes[element] = (startDateMJD[1]+endTimeDelta.days) + long((datetime.min+endTimeDelta).time().hour)/24.0 + \
        long((datetime.min+endTimeDelta).time().minute)/60.0/24.0 + \
        long((datetime.min+endTimeDelta).time().second)/60.0/60.0/24.0+ \
        long((datetime.min+endTimeDelta).time().microsecond)/60.0/60.0/24.0/(1.0e6)
    # The start and end times of the interval are converted to strings in a format ready to be
    # passed to clean, ie. 'hh:mm:ss~hh:mm:ss'.
    if float(endTimeDelta.days) == 0 and float(startTimeDelta.days)==0:
    	timeIntervals[element] = str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +\
        '/' + str(time.isoformat((datetime.min+startTimeDelta).time())) + '~' +str(date_conv[0]) +'/'+\
         str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+endTimeDelta).time()))
    elif float(endTimeDelta.days) != 0 and float(startTimeDelta.days)==0:
    	timeIntervals[element] = str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +\
        '/' + str(time.isoformat((datetime.min+startTimeDelta).time())) + '~' +str(date_conv2[0]) +'/'+\
         str(month_2) +'/'+str(str(date_conv2[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+endTimeDelta).time()))
    elif float(endTimeDelta.days) != 0 and float(startTimeDelta.days)!=0:
    	timeIntervals[element] = str(date_conv2[0]) +'/'+ str(month_2) +'/'+str(str(date_conv2[2]).zfill(2)) +\
        '/' + str(time.isoformat((datetime.min+startTimeDelta).time())) + '~' +str(date_conv2[0]) +'/'+\
         str(month_2) +'/'+str(str(date_conv2[2]).zfill(2)) +'/' + str(time.isoformat((datetime.min+endTimeDelta).time()))
    else:
        raise Exception('Something went wrong finding time intervals,'
                        'please review input')


# If the remainder of observation time is greater than this # of minutes it should be used, and is appended
# to timeIntervals and mjdTimes.
if remainder >= rem_int:
    if float(startDateMJD[1]) != float(endDateMJD[1]):
    	timeIntervals = timeIntervals + [str(date_conv2[0])+'/'+str(month_2)+'/'+str(str(date_conv2[2]).zfill(2)) +\
        '/'+str(time.isoformat((datetime.min+endTimeDelta).time()))+'~'+str(date_conv2[0])+'/'+str(month_2)+'/'\
        +str(str(date_conv2[2]).zfill(2)) +'/'+str(time.isoformat(endTimeHMS))]
    	mjdTimes = mjdTimes + [endDateMJD[1] + long(endTime.hour)/24.0 + long(endTime.minute)/24.0/60.0+ \
        long(endTime.second)/24.0/60.0/60.0+long(endTime.microsecond)/24.0/60.0/60.0/(1.0e6)]
    else:
    	timeIntervals = timeIntervals + [str(date_conv[0]) +'/'+ str(month_0) +'/'+str(str(date_conv[2]).zfill(2)) +\
        '/' +str(time.isoformat((datetime.min+endTimeDelta).time()))+'~'+str(date_conv[0]) +'/'+ str(month_0) +\
        '/'+str(str(date_conv[2]).zfill(2)) +'/' +str(time.isoformat(endTimeHMS))]
    	mjdTimes = mjdTimes + [startDateMJD[1] + long(endTime.hour)/24.0 + long(endTime.minute)/24.0/60.0+ \
        long(endTime.second)/24.0/60.0/60.0+long(endTime.microsecond)/24.0/60.0/60.0/(1.0e6)]

# The results are printed for the user- uncomment if you like
'''print '\nThe observation will be divided into the following intervals: '
print timeIntervals
print '\nmjdTimes'
print mjdTimes'''
##################################

##################################
# CLEAN + Fit each chunk
##################################

#initialize lists
fluxError_real=[]
uv_fitval=[]
uv_fiterr=[]
#make copy of mjdTimes and timeIntervals for uv fitting
mjdTimes_uv=mjdTimes
timeIntervals_uv=timeIntervals
#count failed intervals
counter_fail=0

if runClean == "T":
	print 'CLEAN is starting-->'
	for interval, time, interval_uv, time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
		print 'CLEANing interval: ', interval
		intervalString=interval.replace(':', '.').replace('/','_')
		if outlierFile == '':
			clean(vis=visibility, imagename=outputPath+label+intervalString, timerange=interval,mask=maskPath, selectdata=True, field='', mode='mfs', imsize=imageSize, cell=cellSize, weighting=weighting,robust=robust,spw=spw_choice, nterms=taylorTerms, niter=numberIters, gain=0.1, threshold=thre,interactive=False)
		else:
			clean(vis=visibility, imagename=outputPath+label+intervalString,mask=maskPath, selectdata=True,timerange=interval, field='', mode='mfs', imsize=imageSize, cell=cellSize, weighting=weighting,robust=robust,spw=spw_choice, nterms=taylorTerms, niter=numberIters, gain=0.1, threshold=thre,interactive=False,outlierfile=outlierFile)
		if taylorTerms == 1:
			imSuffix = '.image'
		else:
			imSuffix = '.image.tt0'
		# For some intervals the CLEAN may have failed, in which case the image file for that interval will not exist.
		# An if statement is used to pass over these intervals, avoiding runtime errors.
		if os.path.exists(outputPath+label+intervalString+imSuffix):# and imstat_conv_check['max'][0] <= thre_num:
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
			peakPosX = str(peakPosValue[0])
			peakPosY = str(peakPosValue[1])
			if np.isinf([float(peak),float(peakPosX),float(peakPosY),float(beamMajor.strip('arcsec')),float(beamMinor.strip('arcsec')),float(beamPA.strip('deg'))]).any()==True:
				print '\n Fitting Failed. No fitting done.'
				counter_fail=counter_fail+1
				timeIntervals.remove(interval)
				mjdTimes.remove(time)
				os.system('sudo rm -rf '+outputPath+label+intervalString+'.*')
			else:
				# Save parameters in a temp file which will be passed to imfit. The correct layout is:
				# peak intensity, peak x-pixel value, peak y-pixel value, major axis, minor axis, position angle, fixed
				# the fixed parameter can contain any of the following:
				#'f' (peak intensity), 'x' (peak x position), 'y' (peak y position), 'a' (major axis), 'b' (minor axis), 'p' (position angle)
				#tempFile = tempfile.NamedTemporaryFile()
				tempFile = open('tempfile.txt','w')
				mystring = str(peak+', '+peakPosX+', '+peakPosY+', '+beamMajor+', '+beamMinor+', '+beamPA+',abp')
				tempFile.write(mystring)
				tempFile.close()
			# Run imfit using the above beam parameters.
				if fix_pos== 'F':
					imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox,logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text',estimates=tempFile.name, append=False, overwrite = True)
				elif fix_pos == 'T':
					if do_monte == 'T':
						samp_px=np.random.normal(0,1,nsim)
						samp_py=np.random.normal(0,1,nsim)
						#samp_bmaj=np.random.normal(0,1,nsim)
						#samp_bmin=np.random.normal(0,1,nsim)
						#samp_pos=np.random.normal(0,1,nsim)
						for i in range(0,len(samp_px)):
							peak_x1=(samp_px[i]*peak_x[1])+peak_x[0]
							peak_y1=(samp_py[i]*peak_y[1])+peak_y[0]
							b_maj1=b_maj[0]#(samp_bmaj[i]*b_maj[1])+b_maj[0]
							b_min1=b_min[0]#(samp_bmin[i]*b_min[1])+b_min[0]
							pos_ang1=pos_ang[0]#(samp_pos[i]*pos_ang[1])+pos_ang[0]
							tempFile2 = open('tempfile2.txt','w')
							mystring2 = str(peak+', '+str(peak_x1)+', '+str(peak_y1)+', '+str(b_maj1)+'arcsec, '+str(b_min1)+'arcsec, '+str(pos_ang1)+'deg, '+par_fix)
							tempFile2.write(mystring2)
							tempFile2.close()
							imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox,logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text',estimates=tempFile2.name,append=False, overwrite = True)
					elif do_monte =='F':
						tempFile2 = open('tempfile2.txt','w')
						mystring2 = str(peak+', '+str(peak_x[0])+', '+str(peak_y[0])+', '+str(b_maj[0])+'arcsec, '+str(b_min[0])+'arcsec, '+str(pos_ang[0])+'deg,'+par_fix)
						tempFile2.write(mystring2)
						tempFile2.close()
						imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox,logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text',estimates=tempFile2.name,append=False, overwrite = True)
					else:
						raise Exception('Please specify whether you wish to perform a Monte Carlo fit (T) or not(F)')
				result_box1 = imstat(imagename=outputPath+label+intervalString+imSuffix,region='annulus['+cen_annulus+','+cen_radius+']')
				fluxError_real.append(result_box1['rms'][0])
				if cutout== 'T':
					immath(imagename=outputPath+label+intervalString+imSuffix,mode='evalexpr',expr='IM0',box=cut_reg,outfile=outputPath+label+intervalString+'_temp'+imSuffix)
					immath(imagename=outputPath+label+intervalString+imSuffix,mode='evalexpr',expr='IM0',region='annulus['+cen_annulus+','+cen_radius+']',outfile=outputPath+label+intervalString+'_rms'+imSuffix)
					os.system('sudo rm -rf '+outputPath+label+intervalString+'.*')
					os.system('sudo mv '+outputPath+label+intervalString+'_temp'+imSuffix+' '+outputPath+label+intervalString+imSuffix)
		else:
			print '\nCLEAN failed on interval ' + interval + '.'
			counter_fail=counter_fail+1
			# The corresponding time intervals must be removed from timeIntervals to avoid runtime errors further on.
			timeIntervals.remove(interval)
			mjdTimes.remove(time)
			os.system('sudo rm -rf '+outputPath+label+intervalString+'.*')
		#uvfitting if requested
		if uv_fit =='T':
			if np.where(np.array(timeIntervals)==interval)[0][0]==0:
				combuv=visibility_uv.strip('.ms')
				mstransform(vis=visibility_uv,outputvis=combuv+'_mstransform.ms',combinespws=True,spw='',datacolumn='data')
				fit=uvm.uvmultifit(vis=combuv+'_mstransform.ms', MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv],spw=spw_choice, column = "data", uniform=False, model=[comp_uv],stokes = stokes_param,var=[uv_var], p_ini=[float(src_uv_init.split(',')[0]),float(src_uv_init.split(',')[1]),float(src_uv_init.split(',')[2])],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.txt',OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
			else:
				fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
				fit.fit(redo_fixed=False,reinit_model=False)
			uv_fitval.append(fit.result['Parameters'][2])
			uv_fiterr.append(fit.result['Uncertainties'][2])
elif runClean == "F":
	print 'Fitting is starting-->'
	counter_fail=0
	for interval, time, interval_uv,time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
		print 'Fitting interval: ', interval
		intervalString=interval.replace(':', '.').replace('/','_')
		if outlierFile == '':
			if taylorTerms == 1:
				imSuffix = '.image'
			else:
				imSuffix = '.image.tt0'
		else:
			if taylorTerms == 1:
				imSuffix = '.image'
			else:
				imSuffix = '.image.tt0'
		if os.path.exists(outputPath+label+intervalString+imSuffix):# and imstat_conv_check['max'][0] <= thre_num:
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
			peakPosX = str(peakPosValue[0])
			peakPosY = str(peakPosValue[1])
			if np.isinf([float(peak),float(peakPosX),float(peakPosY),float(beamMajor.strip('arcsec')),float(beamMinor.strip('arcsec')),float(beamPA.strip('deg'))]).any()==True:
				print '\n Fitting Failed. No fitting done.'
				counter_fail=counter_fail+1
				timeIntervals.remove(interval)
				mjdTimes.remove(time)
			else:
				tempFile = open('tempfile.txt','w')
				mystring = str(peak+', '+peakPosX+', '+peakPosY+', '+beamMajor+', '+beamMinor+', '+beamPA+',abp')
				tempFile.write(mystring)
				tempFile.close()
				if fix_pos=='F':
					imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox,logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text',estimates=tempFile.name, append=False, overwrite = True)
				else:
					if do_monte == 'T':
						samp_px=np.random.normal(0,1,nsim)
						samp_py=np.random.normal(0,1,nsim)
						#samp_bmaj=np.random.normal(0,1,nsim)
						#samp_bmin=np.random.normal(0,1,nsim)
						#samp_pos=np.random.normal(0,1,nsim)
						for i in range(0,len(samp_px)):
							peak_x1=(samp_px[i]*peak_x[1])+peak_x[0]
							peak_y1=(samp_py[i]*peak_y[1])+peak_y[0]
							b_maj1=b_maj[0]
							b_min1=b_min[0]
							pos_ang1=pos_ang[0]
							tempFile2 = open('tempfile2.txt','w')
							mystring2 = str(peak+', '+str(peak_x1)+', '+str(peak_y1)+', '+str(b_maj1)+'arcsec, '+str(b_min1)+'arcsec, '+str(pos_ang1)+'deg, '+par_fix)
							tempFile2.write(mystring2)
							tempFile2.close()
							imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox,logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text',estimates=tempFile2.name,append=False, overwrite = True)
					elif do_monte =='F':
						tempFile2 = open('tempfile2.txt','w')
						mystring2 = str(peak+', '+str(peak_x[0])+', '+str(peak_y[0])+', '+str(b_maj[0])+'arcsec, '+str(b_min[0])+'arcsec, '+str(pos_ang[0])+'deg,'+par_fix)
						tempFile2.write(mystring2)
						tempFile2.close()
						imfit(imagename=outputPath+label+intervalString+imSuffix, box=targetBox,logfile=outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text',estimates=tempFile2.name,append=False, overwrite = True)
					else:
						raise Exception('Please specify whether you wish to perform a Monte Carlo fit (T) or not(F)')
				result_box1=imstat(imagename=outputPath+label+intervalString+imSuffix,region='annulus['+cen_annulus+','+cen_radius+']')
				fluxError_real.append(result_box1['rms'][0])
		else:
			print '\nCLEAN failed on interval ' + interval + '. No fitting done'
			# The corresponding time intervals must be removed from timeIntervals to avoid runtime errors further on.
			counter_fail=counter_fail+1
			timeIntervals.remove(interval)
			mjdTimes.remove(time)
		if uv_fit =='T':
			if np.where(np.array(timeIntervals)==interval)[0][0]==0:
				combuv=visibility_uv.strip('.ms')
				mstransform(vis=visibility_uv,outputvis=combuv+'_mstransform.ms',combinespws=True,spw='',datacolumn='data')
				fit=uvm.uvmultifit(vis=combuv+'_mstransform.ms',MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv],spw=spw_choice, column = "data", uniform=False, model=[comp_uv],stokes = stokes_param,var=[uv_var], p_ini=[float(src_uv_init.split(',')[0]),float(src_uv_init.split(',')[1]),float(src_uv_init.split(',')[2])],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.txt',OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
			else:
				fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
				fit.fit(redo_fixed=False,reinit_model=False)
			uv_fitval.append(fit.result['Parameters'][2])
			uv_fiterr.append(fit.result['Uncertainties'][2])
elif runClean == 'U':
	print 'You selected only UV plane analysis. No image plane analysis will be done.'
	for interval, time, interval_uv, time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
		print 'UV analysis on interval: ', interval
		intervalString=interval.replace(':', '.').replace('/','_')
		if uv_fit =='T':
			if np.where(np.array(timeIntervals)==interval)[0][0]==0:
				combuv=visibility_uv.strip('.ms')
				mstransform(vis=visibility_uv,outputvis=combuv+'_mstransform.ms',combinespws=True,spw='',datacolumn='data')
				fit=uvm.uvmultifit(vis=combuv+'_mstransform.ms', MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv],spw=spw_choice, column = "data", uniform=False, model=[comp_uv],stokes = stokes_param,var=[uv_var], p_ini=[float(src_uv_init.split(',')[0]),float(src_uv_init.split(',')[1]),float(src_uv_init.split(',')[2])],outfile = outputPath+'uvfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.txt',OneFitPerChannel=False ,cov_return=False, finetune=False, method="levenberg")
			else:
				fit.MJDrange=[time_uv-(intervalSizeH/24.+intervalSizeM/(24.*60.)+intervalSizeS/(24.*60.*60.)),time_uv]
				fit.fit(redo_fixed=False,reinit_model=False)
			uv_fitval.append(fit.result['Parameters'][2])
			uv_fiterr.append(fit.result['Uncertainties'][2])

uv_fitval_arr=np.array(uv_fitval)
uv_fiterr_arr=np.array(uv_fiterr)
##################################

##################################
# Creating output data arrays
##################################
print 'Creating output flux arrays...\n'
fluxDensity = []#integ
fluxDensity2 =[]#peak
fluxDensity3 =[]#UV plane
fluxError = []#integ
fluxError2 =[]#peak
fluxError3 =[]#UV plane
suffix = []#integ
suffix2 = []#peak
suffix3 = []#UV plane
timerange = []

# Create flux density, flux density error, suffix and timerange vectors. In each case get values from imfit text files.
# The exact layout of the imfit outputs are somewhat unpredictable. Therefore find() is used here to reliably locate
# the relavent information.
for interval, time,interval_uv,time_uv in zip(timeIntervals, mjdTimes,timeIntervals_uv, mjdTimes_uv):
	flux_uv=[]
	fluxerr_uv=[]
	if runClean != 'U':
		if do_monte =='T':
			FD_list=[]
			FDE_list=[]
			FD2_list=[]
			FDE2_list=[]
			fluxerr_uv=[]
			for i in range(0,len(samp_px)):
				intervalString=interval.replace(':', '.').replace('/','_')
				if os.path.exists(outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text'):
					imfitFile = open(outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+str(i)+'.text', 'r')
					imfitText = imfitFile.read()
					imfitFile.close()
					end = imfitText[-19:-1]
					if end != '*** FIT FAILED ***':
						if integ_fit == 'T':
							startWord = '--- Integrated:'
							endWord = '--- Peak:'
							startPos = imfitText.find(startWord)
							endPos = imfitText.find(endWord)
							fluxString = imfitText[startPos:endPos]
							plusMinusPos = fluxString.find('+/-')
							unitsPos = fluxString.find('Jy')
							plot_label_unit=''
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
							else:
								fluxD=fluxD*[1.0]
								fluxE=fluxE*[1.0]
							FD_list.append(fluxD)
							FDE_list.append(fluxE)
						elif integ_fit == 'F':
							startWord = '--- Peak:'
							endWord = '--- Polarization:'
							startPos = imfitText.find(startWord)
							endPos = imfitText.find(endWord)
							fluxString = imfitText[startPos:endPos]
							plot_label_unit='/beam'
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
							else:
								fluxD=fluxD*1.0
								fluxE=fluxE*1.0
							FD_list.append(fluxD)
							FDE_list.append(fluxE)
						elif integ_fit == 'B':
							startWord = '--- Integrated:'
							startWord2 = '--- Peak:'
							endWord = '--- Peak:'
							endWord2 = '--- Polarization:'
							startPos = imfitText.find(startWord)
							startPos2 = imfitText.find(startWord2)
							endPos = imfitText.find(endWord)
							endPos2 = imfitText.find(endWord2)
							fluxString = imfitText[startPos:endPos]
							fluxString2 = imfitText[startPos2:endPos2]
							plusMinusPos = fluxString.find('+/-')
							plusMinusPos2 = fluxString2.find('+/-')
							unitsPos = fluxString.find('Jy')
							unitsPos2 = fluxString2.find('Jy/b')
							plot_label_unit=''
							plot_label_unit2='/beam'
							fluxD1=float(fluxString[16:plusMinusPos])
							fluxE1=float(fluxString[plusMinusPos+3:unitsPos-1])
							fluxD2=float(fluxString2[10:plusMinusPos2])
							fluxE2=float(fluxString2[plusMinusPos2+3:unitsPos2-1])
							suff = str(fluxString[unitsPos-1])
							suff2 = str(fluxString2[unitsPos2-1])
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
							else:
								fluxD1=fluxD1*1.0
								fluxE1=fluxE1*1.0
							if suff2 == 'u':
								fluxD2=fluxD2*1.0e-6
								fluxE2=fluxE2*1.0e-6
							elif suff2 == 'm':
								fluxD2=fluxD2*1.0e-3
								fluxE2=fluxE2*1.0e-3
							elif suff2 == 'n':
								fluxD2=fluxD2*1.0e-9
								fluxE2=fluxE2*1.0e-9
							else:
								fluxD2=fluxD2*1.0
								fluxE2=fluxE2*1.0

							FD_list.append(fluxD1)
							FDE_list.append(fluxE1)
							FD2_list.append(fluxD2)
							FDE2_list.append(fluxE2)
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
			timerange = timerange + [interval]
			suffix = suffix + ['m']
		elif do_monte == 'F':
			intervalString=interval.replace(':', '.').replace('/','_')
			if os.path.exists(outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text'):
				imfitFile = open(outputPath+'imfit_'+target+'_'+refFrequency+'_'+obsDate+'_'+intervalString+'.text', 'r')
				imfitText = imfitFile.read()
				imfitFile.close()
				end = imfitText[-19:-1]
				if end != '*** FIT FAILED ***':
					if integ_fit == 'T':
						startWord = '--- Integrated:'
						endWord = '--- Peak:'
						startPos = imfitText.find(startWord)
						endPos = imfitText.find(endWord)
						fluxString = imfitText[startPos:endPos]
						plusMinusPos = fluxString.find('+/-')
						unitsPos = fluxString.find('Jy')
						plot_label_unit=''
						fluxDensity = fluxDensity + [float(fluxString[16:plusMinusPos])]
						fluxError = fluxError + [float(fluxString[plusMinusPos+3:unitsPos-1])]
						suffix = suffix + [str(fluxString[unitsPos-1])]
					elif integ_fit == 'F':
						startWord = '--- Peak:'
						endWord = '--- Polarization:'
						startPos = imfitText.find(startWord)
						endPos = imfitText.find(endWord)
						fluxString = imfitText[startPos:endPos]
						plot_label_unit='/beam'
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
						plot_label_unit=''
						fluxDensity= fluxDensity +[float(fluxString[16:plusMinusPos])]
						fluxError=fluxError +[float(fluxString[plusMinusPos+3:unitsPos-1])]
						suffix = suffix + [str(fluxString[unitsPos-1])]
						startWord2 = '--- Peak:'
						endWord2 = '--- Polarization:'
						startPos2 = imfitText.find(startWord2)
						endPos2 = imfitText.find(endWord2)
						fluxString2 = imfitText[startPos2:endPos2]
						plot_label_unit2='/beam'
						plusMinusPos2 = fluxString2.find('+/-')
						unitsPos2 = fluxString2.find('Jy/b')
						fluxDensity2=fluxDensity2+[float(fluxString2[10:plusMinusPos2])]
						fluxError2=fluxError2+[float(fluxString2[plusMinusPos2+3:unitsPos2-1])]
						suffix2 = suffix2 + [str(fluxString2[unitsPos2-1])]
					timerange = timerange + [interval]
				else:
					print '\nFit failed on interval ' + interval
					fluxError_real.remove(fluxError_real[timeIntervals.index(interval)])
					timeIntervals.remove(interval)
					mjdTimes.remove(time)
			else:
				fluxError_real.remove(fluxError_real[timeIntervals.index(interval)])
				timeIntervals.remove(interval)
				mjdTimes.remove(time)
		else:
			raise Exception('Please specify whether you wish to perform a Monte Carlo fit (T) or not(F)')
		if uv_fit == 'T':
			flux_uv1=uv_fitval_arr[np.where(np.array(timeIntervals)==interval)[0]][0]
			fluxerr_uv1=uv_fiterr_arr[np.where(np.array(timeIntervals)==interval)[0]][0]
			fluxDensity3 = fluxDensity3 + [flux_uv1]
			fluxError3 = fluxError3 + [fluxerr_uv1]
			suffix3 = suffix3 + [' ']
	elif runClean=='U':
		if uv_fit == 'T':
			flux_uv1=uv_fitval_arr[np.where(np.array(timeIntervals)==interval)[0]][0]
			fluxerr_uv1=uv_fiterr_arr[np.where(np.array(timeIntervals)==interval)[0]][0]
			fluxDensity3 = fluxDensity3 + [flux_uv1]
			fluxError3 = fluxError3 + [fluxerr_uv1]
			suffix3 = suffix3 + [' ']
##################################


##################################
#Print raw results to terminal- uncomment if you like
##################################
'''print 'Printing raw results for user...\n'
if runClean != 'U':
	if integ_fit == 'B':
		if len(fluxDensity)==0 or fluxDensity2==0:
			raise Exception('CLEAN failed on all intervals. Please review input.')
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
		if len(fluxDensity)==0:
			raise Exception('CLEAN failed on all intervals. Please review input.')
		print '\nFlux desities extracted from imfit:'
		print fluxDensity
		print '\nFlux errors extracted from imfit:'
		print fluxError
		print '\nUnits extracted from imfit:'
		print suffix
		print '\nRealistic Error from regions (Jy):'
		print fluxError_real

	if uv_fit == 'T':
		if len(fluxDensity3)==0:
			raise Exception('UV fitting failed on all intervals. Please review input.')
		print '\nFlux desities extracted from uvmodelfit:'
		print fluxDensity3
		print '\nFlux errors extracted from uvmodelfit:'
		print fluxError3
		print '\nUnits extracted from uvmodelfit:'
		print suffix3
else:
	if len(fluxDensity3)==0:
		raise Exception('UV fitting failed on all intervals. Please review input.')
	print '\nFlux desities extracted from uvmodelfit:'
	print fluxDensity3
	print '\nFlux errors extracted from uvmodelfit:'
	print fluxError3
	print '\nUnits extracted from uvmodelfit:'
	print suffix3'''
##################################


##################################
#Convert to consistent units
##################################
print 'Checking flux units...\n'
# The flux values will be plotted in a user specified unit, however the imfit outputs will not necessarily be in these
# units. To get around this the units for each value are recorded in the suffix list created in the above
# for loop. The following loop will use this list to create a multiplier list, which will be used to
# convert all values to Jy. They can then all be scaled to be ready for plotting.
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
    else:
        multiplier = multiplier + [1.0]
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
        else:
            multiplier2 = multiplier2 + [1.0]
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
        else:
            multiplier3 = multiplier3 + [1.0]

# Multiplier() should now contain a list of multiplication factors ready to be applied to the flux data. The
#end result is in Jy
if runClean=='U':
	length= range(len(fluxDensity3))
	for n in length:
		fluxDensity3[n] = fluxDensity3[n] * multiplier3[n]
		fluxError3[n] = fluxError3[n] * multiplier3[n]
else:
	length = range(len(fluxDensity))
	for n in length:
		fluxDensity[n] = fluxDensity[n] * multiplier[n]
		fluxError[n] = fluxError[n] * multiplier[n]
		if integ_fit == 'B':
		    fluxDensity2[n] = fluxDensity2[n] * multiplier2[n]
		    fluxError2[n] = fluxError2[n] * multiplier2[n]
		if uv_fit == 'T':
		    fluxDensity3[n] = fluxDensity3[n] * multiplier3[n]
		    fluxError3[n] = fluxError3[n] * multiplier3[n]

# The flux values will be plotted in a user specified unit, so they are scaled here,
#along with error values.
if runClean=='U':
	length= range(len(fluxDensity3))
	for k in length:
		fluxDensity3[k] = fluxDensity3[k]*lc_scale_factor
		fluxError3[k] = fluxError3[k]*lc_scale_factor
else:
	length= range(len(fluxDensity))
	for k in length:
		fluxDensity[k] = fluxDensity[k]*lc_scale_factor
		fluxError[k] = fluxError[k]*lc_scale_factor
		fluxError_real[k]=fluxError_real[k]*lc_scale_factor
		if integ_fit == 'B':
		    fluxDensity2[k] = fluxDensity2[k]*lc_scale_factor
		    fluxError2[k] = fluxError2[k]*lc_scale_factor
		if uv_fit == 'T':
		    fluxDensity3[k] = fluxDensity3[k]*lc_scale_factor
		    fluxError3[k] = fluxError3[k]*lc_scale_factor
##################################

if runClean=='U':
	if len(fluxDensity3)==0:
		print 'Fitting failed in all timebins. Please check input.'
		failed='y'
	else:
		failed='n'
else:
	if len(fluxDensity)==0:
		print 'Fitting failed in all timebins. Please check input.'
		failed='y'
	else:
		failed='n'
##################################
#Write results to data file
##################################
if failed=='n':
	print 'Writing light curve data file...\n'
	data = open(dataPath, 'w')
	if runClean != 'U':
	    #data.write('{0} {1}\n'.format('#The number of timebins where CLEAN failed/did not converge was',counter_fail))
	    if integ_fit == 'B':
	    	data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|peak flux|peak imfit err|rms error'))
	    elif integ_fit == 'B' and uv_fit == 'T':
	    	data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|peak flux|peak imfit err|uv flux|uverr|rms error'))
	    elif integ_fit=='T':
	    	data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|rms error'))
	    elif integ_fit=='T' and uv_fit=='T':
	    	data.write('{0}\n'.format('#MJD|integ flux|integ imfit err|uv flux|uverr|rms error'))
	    elif integ_fit=='F':
	    	data.write('{0}\n'.format('#MJD|peak flux|peak imfit err|rms error'))
	    elif integ_fit=='F' and uv_fit=='T':
	    	data.write('{0}\n'.format('#MJD|peak flux|peak imfit err|uv flux|uverr|rms error'))
	    for i in range(0,len(fluxDensity)):
	        if integ_fit == 'B':
	            data.write('{0:.12f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity2[i],\
	                fluxError2[i],fluxError_real[i]))
	            if uv_fit == 'T':
	                data.write('{0:.12f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity2[i],\
	                    fluxError2[i],fluxDensity3[i],fluxError3[i],fluxError_real[i]))
	        elif integ_fit =='T':
	            data.write('{0:.12f} {1:.3f} {2:.3f} {3:.3f}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxError_real[i]))
	            if uv_fit=='T':
	                data.write('{0:.12f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity3[i],\
	                    fluxError3[i],fluxError_real[i]))
	        elif integ_fit =='F':
	            data.write('{0:.12f} {1:.3f} {2:.3f} {3:.3f}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxError_real[i]))
	            if uv_fit=='T':
	                data.write('{0:.12f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}\n'.format(mjdTimes[i],fluxDensity[i],fluxError[i],fluxDensity3[i],\
	                    fluxError3[i],fluxError_real[i]))
	else:
	    data.write('{0}\n'.format('#MJD|uv flux|uverr'))
	    for i in range(0,len(fluxDensity3)):
	        data.write('{0:.12f} {1:.3f} {2:.3f}\n'.format(mjdTimes[i],fluxDensity3[i],fluxError3[i]))
	data.close()
	print dataPath+' is saved.'
	##################################
	##################################
	#Plot Lightcurves
	##################################
	print 'Plotting Light curves...\n'
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
	if runClean != 'U':
	    if integ_fit == 'B':
	    	fig1=pp.figure()
	    	pp.errorbar(Elapsed, fluxDensity, yerr=fluxError_real, fmt='ro',)
	    	pp.xlabel('Time since start of observation (mins)')
	    	y_label_name='Flux Density (mJy'+plot_label_unit+')'
	    	pp.ylabel(y_label_name)
	    	pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	    	pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	    	savestring = os.path.join(path_dir,
	                                  'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+
	                                  str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_integ.pdf')
	    	pp.savefig(savestring)
	    	print savestring, ' is saved'
	    	fig2=pp.figure()
	    	pp.errorbar(Elapsed, fluxDensity2, yerr=fluxError_real, fmt='ro',)
	    	pp.xlabel('Time since start of observation (mins)')
	    	y_label_name2='Flux Density (mJy'+plot_label_unit2+')'
	    	pp.ylabel(y_label_name2)
	    	pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	    	pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	    	savestring2 = os.path.join(path_dir, 'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+ str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_peak.pdf')
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
	    		savestring = os.path.join(path_dir, 'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_check_lc_uv.pdf')
	    		pp.savefig(savestring)
	    		print savestring, ' is saved'

	    else:
	    	pp.errorbar(Elapsed, fluxDensity, yerr=fluxError_real, fmt='ro',)
	    	pp.xlabel('Time since start of observation (mins)')
	    	y_label_name='Flux Density (mJy'+plot_label_unit+')'
	    	pp.ylabel(y_label_name)
	    	pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	    	pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	    	savestring = os.path.join(path_dir, 'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'.pdf')
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
	    		savestring = os.path.join(path_dir, 'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_uv.pdf')
	    		pp.savefig(savestring)
	    		print savestring, ' is saved'
	else:
	    fig4=pp.figure()
	    pp.errorbar(Elapsed, fluxDensity3, yerr=fluxError3, fmt='ro',)
	    pp.xlabel('Time since start of observation (mins)')
	    y_label_name='Flux Density (mJy/beam)'
	    pp.ylabel(y_label_name)
	    pp.title('Flux Density vs Time. '+target+' '+refFrequency)
	    pp.xlim(0, Elapsed[len(Elapsed)-1]+intervalSizeS/(60.0))
	    savestring = os.path.join(path_dir, 'data_products/'+target+lab+str(intervalSizeH)+'hour_'+str(intervalSizeM)+'min_'+\
	        str(intervalSizeS)+'sec_'+refFrequency+'_'+obsDate+'_uv.pdf')
	    pp.savefig(savestring)
	    print savestring, ' is saved'
	##################################
	##################################
	#Basic Variability Tests
	##################################
	if var_anal=='T':
		print 'Performing Variability Analysis...\n'
		if runClean =='U':
			fluxvar=np.array(fluxDensity3)
			fluxerrvar=np.array(fluxError3)
		else:
			fluxerrvar=np.array(fluxError_real)
			if integ_fit=='T':
				fluxvar=np.array(fluxDensity)
			elif integ_fit=='F':
				fluxvar=np.array(fluxDensity2)
			else:
				fluxvar=np.array(fluxDensity2)
		var_file=open(dataPathVar,'w')
		print 'Performing Variiability Tests'
		chi_tot,dof,null,wm,wmerr,ex_var,ex_var_error,frac_rms,frac_rms_error=var_analysis(fluxvar,fluxerrvar)
		if power_spec=='T':
			print 'Creating Power Spectrum'
			if float(intervalSizeDelta.seconds)==0.0:
				sig95,sig99=lomb_scargle(mjdTimes,fluxvar,fluxerrvar,float(intervalSizeDelta.microseconds)/1e6,labelP,'lin')
			else:
				sig95,sig99=lomb_scargle(mjdTimes,fluxvar,fluxerrvar,float(intervalSizeDelta.seconds),labelP,'lin')
			print labelP+' is saved.'
		var_file.write('{0} {1} {2}\n'.format('Weighted Mean/Error',wm,wmerr))
		var_file.write('{0} {1} {2}\n'.format('Chi2 with weighted mean/dof',chi_tot,dof))
		var_file.write('{0} {1} {2}\n'.format('Excess Variance/Error',ex_var,ex_var_error))
		var_file.write('{0} {1} {2}\n'.format('Fractional RMS/Error',frac_rms,frac_rms_error))
		if power_spec=='T':
			var_file.write('{0} {1} {2}\n'.format('95% and 99% Significance Levels for Periodogram',sig95,sig99))
		var_file.close()
		print dataPathVar+ 'is saved.'
	##################################

##################################
#cleaning up
##################################
print 'Cleaning up...\n'
#remove temp files and .last/.log files created by CASA/ipython
os.system('sudo rm -rf *.last')
os.system('sudo rm -rf casa*.log')
os.system('sudo rm -rf tempfile.txt tempfile2.txt')
print '*********************************************************'
print 'Script finished. Please inspect resulting data products'
print '*********************************************************'
##################################