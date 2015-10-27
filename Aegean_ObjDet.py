############################################################################################################
#Object Detection script-->Using AEGEAN algorithm (https://github.com/PaulHancock/Aegean)
#Input: Cleaned FITS image of whole data set
#Output: Data file of the parameters of objects found in image
#NOTE: Run in normal python, DO NOT RUN WITHIN CASA, it won't work.             
#############################################################################################################
#Written by A. Tetarenko--> 10/2015
#############################################################################################################
#Import modules                                         
#Make sure aegean tree of directories in path so you can import it
import aegean
from AegeanTools.catalogs import save_catalog
import numpy as np
import multiprocessing

#set path to where output is to be stored-->need to set up file system so have data and data_products directory
#in this path
path_dir='/home/ubuntu/'

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
outputPath = path_dir+'data_products/images_'+target+'_'+refFrequency+'_'+str(intervalSizeH)+'hours'+str(intervalSizeM)+'min'+str(intervalSizeS)+'sec/'

fits_file=outputPath+label+'whole_dataset.fits'
out_file0=outputPath+label+'whole_dataset_aegean.txt'
tab_file=outputPath+label+'whole_dataset_objdet.tab'
seed=10
flood=4
tele='VLA'

###########################################################
#END OF USER INPUT SECTION
###########################################################


#use Aegean object detection algorithm--> import as module
ra_list=[]
dec_list=[]
maj_list=[]
min_list=[]
pos_list=[]
src_list=[]
sources=[]
lat=aegean.scope2lat(tele)
if lat==None:
    if tele=='SMA':
        lat=19.8243
    elif tele=='NOEMA':
        lat=44.6339
    else:
        lat_string=raw_input('Enter latitude of telescope-->')
        lat=float(lat_string))
out_file= open(out_file0, 'w')
tables=tab_file
print 'Running Aegean Object Detection -->'
detections=aegean.find_sources_in_image(fits_file, outfile=out_file, hdu_index=0,
                                           rms=None,
                                           max_summits=None,
                                           innerclip=seed,
                                           outerclip=flood, cores=multiprocessing.cpu_count(), rmsin=None,
                                           bkgin=None, beam=None,
                                           doislandflux=False,
                                           nonegative=not False, nopositive=False,
                                           mask=None, lat=lat, imgpsf=None)
out_file.close()
if len(detections) == 0:
    print 'No sources detected'
sources.extend(detections)
if len(sources) > 0:
    save_catalog(tables, sources)
