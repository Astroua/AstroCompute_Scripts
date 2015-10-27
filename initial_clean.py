############################################################################################################
#CASA script #1-->initial CLEAN of FUll data set
#Input: Calibrated and Split MS
#Output: Cleaned and deconvolved CASA image (.image,.flux,.psf,.residual,.mask,.model) & fits image
#Note: This script is theoretically compatible with any data that can be imported 
#into CASA, but has only been tested on continuum data from the VLA, SMA, and NOEMA
#(import NOEMA data into CASA: http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf).             
#############################################################################################################
#Written by A. Tetarenko--> 10/2015
#############################################################################################################
#Import modules                                         
#
#None needed

#set path to where output is to be stored-->need to set up file system so have data and data_products directory
#in this path
path_dir='/home/ubuntu/'

#Function to find next power of 2^n closest to chosen imsize value to optimize cleaning
def is_power2(num):
	return(num !=0 and ((num & (num-1)) ==0))


###########################################################
#USER INPUT SECTION--> will make parameter file eventually
###########################################################
# Target name.
target = 'V404Cyg'
# Date of observation.
obsDate = '2015jun22'
# Observation frequency. 
refFrequency ='21GHz'
# Label for casa output directories and files.
label = target + '_' + refFrequency + '_' + obsDate + '_'
# Length of time bins (H,M,S); see below if you want manual input (line 288)
intervalSizeH = 0
intervalSizeM = 2
intervalSizeS = 0
#make data_products directory before start script
mkdir_string='sudo mkdir '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
mkdir_perm1='sudo chown ubuntu '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
mkdir_perm2='sudo chmod -R 777 '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
os.system(mkdir_string)
os.system(mkdir_perm1)
os.system(mkdir_perm2)
# Path to directory where all output from this script is saved.
outputPath = path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec/'
# dataPath contains the path and filename in which data file will be saved. 
# This script can be run on several epochs of data from the same observation without changing this path.
# In this case the data file will be appended each time.
dataPath = path_dir+'data_products/datafile_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec.txt'
# Name of visibliity - should include full path if script is not being run from vis location.
visibility = path_dir+'data/swj17_jun22_B_K_k21.ms'

# The clean command (line 425) should be inspected closely to ensure all arguments are appropriate before 
# running script on a new data set.
# The following arguments will be passed to casa's clean, imfit or imstat functions:
imageSize = [6000,6000]
if is_power2(imageSize[0])==False:
	print 'Clean will run faster if image size is 2^n'
	imS=raw_input('Do you want to optimize image size?, y or n? ')
	if imS == 'y':
		imageSize=[int(pow(2, m.ceil(np.log(imageSize[0])/np.log(2)))),int(pow(2, m.ceil(np.log(imageSize[0])/np.log(2))))]
        	print 'imagesize is now set to ', imageSize
    	print 'imagesize remains at ',imageSize, 'due to user request'
numberIters = 5000
cellSize = ['0.02arcsec','0.02arcsec']
taylorTerms = 1
myStokes = 'I'
thre='10mJy'
spw_choice='0~7:5~58'
###########################################################
#END OF USER INPUT SECTION
###########################################################

#CLEAN and CONVERT to FITS for use with object detection script
print 'Cleaning Full Data Set to detect objects-->'
clean(vis=visibility, imagename=outputPath+label+'whole_dataset', field='', mode='mfs', imsize=imageSize, cell=cellSize, weighting='natural',spw=spw_choice, nterms=taylorTerms, niter=numberIters, gain=0.1, threshold=thre, interactive=F)
print 'Converting to fits-->'
exportfits(imagename=outputPath+label+'whole_dataset.image', fitsimage=outputPath+label+'whole_dataset.fits',history=False)
