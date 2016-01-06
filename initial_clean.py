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

import json
import os
import sys

sys.path.append("AstroCompute_Scripts/")
from AstroCompute_Scripts.utils import load_json

#set path to where output is to be stored-->need to set up file system so have data and data_products directory
#in this path
path_dir = sys.argv[-1]
# path_dir='/home/ubuntu/'

param_file = sys.argv[-2]

###########################################################
#USER INPUT SECTION-->read in from parameters file
###########################################################
data_params = load_json(param_file)

# Target name.
target = data_params["target"]
# Date of observation.
obsDate = data_params["obsDate"]
# Observation frequency.
refFrequency = data_params["refFrequency"]
# Label for casa output directories and files.
label = target + '_' + refFrequency + '_' + obsDate + '_'
# Length of time bins (H,M,S)
intervalSizeH = data_params["intervalSizeH"]
intervalSizeM = data_params["intervalSizeM"]
intervalSizeS = data_params["intervalSizeS"]
#make data_products directory before start script
mkdir_string='sudo mkdir '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
mkdir_perm1='sudo chown ubuntu '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
mkdir_perm2='sudo chmod -R 777 '+path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec'
os.system(mkdir_string)
os.system(mkdir_perm1)
os.system(mkdir_perm2)
# Path to directory where all output from this script is saved.
outputPath = \
    os.path.join(path_dir, 'data_products/images_'+target+'_'+refFrequency +
                 '_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min' +
                 str(intervalSizeS)+'sec/')
# dataPath contains the path and filename in which data file will be saved.
# This script can be run on several epochs of data from the same observation without changing this path.
# In this case the data file will be appended each time.
dataPath = path_dir+'data_products/datafile_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec.txt'
# Name of visibliity - should include full path if script is not being run from vis location.
visibility = path_dir+'data/'+ data_params["visibility"]

# The following arguments will be passed to casa's clean
imageSize = [data_params["imageSize"]] * 2
numberIters = data_params["numberIters"]
cellSize = [data_params["cellSize"]]
taylorTerms = data_params["taylorTerms"]
myStokes = data_params["myStokes"]
thre=data_params["thre"]
spw_choice=data_params["spw_choice"]
###########################################################
#END OF USER INPUT SECTION
###########################################################

#CLEAN and CONVERT to FITS for use with object detection script
print 'Cleaning Full Data Set to detect objects-->'
clean(vis=visibility,
      imagename=os.path.join(outputPath, label+'whole_dataset'),
      field='', mode='mfs', imsize=imageSize, cell=cellSize,
      weighting='natural', spw=spw_choice, nterms=taylorTerms,
      niter=numberIters, gain=0.1,
      threshold=thre, interactive=False)
print 'Converting to fits-->'
exportfits(imagename=os.path.join(outputPath, label+'whole_dataset.image'),
           fitsimage=os.path.join(outputPath, label+'whole_dataset.fits'),
           history=False)
