##################################
#Object Detection Script
##################################
'''Uses AEGEAN algorithm (https://github.com/PaulHancock/Aegean) for object detection in a radio image.
INPUT: Cleaned FITS image of whole data set--> [fits_file]
OUTPUT: Data file of the parameters of objects found in image--> [tab_file]
NOTE: Needs to be run with lmfit 0.7.4, otherwise it won't work.
(casa-pip install git+git://github.com/lmfit/lmfit-py.git@0.7.4)

Written by: A. Tetarenko
Last Updated: February 2 2017'''

# Import modules
# Make sure aegean tree of directories in path so you can import it
import sys
sys.path.append('/home/ubuntu/Aegean')
import aegean
from AegeanTools.catalogs import save_catalog
import numpy as np
from multiprocessing import cpu_count
import re
import os
from utils import run_aegean,initial_clean
import warnings
warnings.filterwarnings('ignore')

def objdet(tele,lat,out_file0,fits_file,seed,flood,tab_file,catalog_input_name,cellSize_string):
  sources = []

  # get latitude for telescope
  lat = aegean.scope2lat(tele)
  # adding SMA and NOEMA to list
  if lat == None:
    if tele == 'SMA':
        lat = 19.8243
    elif tele == 'NOEMA':
        lat = 44.6339
    else:
        #lat_string=raw_input('Enter latitude of telescope-->')
        lat = data_params["lat"]

  out_file = open(out_file0, 'w')

  print 'Running Aegean Object Detection -->'
  detections = aegean.find_sources_in_image(fits_file,
                                          outfile=out_file,
                                          hdu_index=0,
                                          rms=None,
                                          max_summits=None,
                                          innerclip=seed,
                                          outerclip=flood,
                                          cores=cpu_count(),
                                          rmsin=None,
                                          bkgin=None, beam=None,
                                          doislandflux=False,
                                          nonegative=not False,
                                          nopositive=False,
                                          mask=None, lat=lat, imgpsf=None)
  out_file.flush()
  out_file.close()
  if len(detections) == 0:
    raise Exception('No sources detected by Aegean. Please check your inputs.')
  sources.extend(detections)
  # write detected source info to file
  if len(sources) > 0:
    save_catalog(catalog_input_name, sources)

  # print detected sources fro user to choose from
  src_l, ra_l, dec_l, maj_l, min_l, pos_l = run_aegean(tab_file, cellSize_string)
  print 'Number of Objects Detected is ', len(src_l)
  print 'Objects Detected-->'
  print 'Object, RA, DEC'
  for i in range(0, len(src_l)):
    print src_l[i], ra_l[i], dec_l[i]
  return(src_l, ra_l, dec_l, maj_l, min_l, pos_l)

if __name__ == "__main__":
  ##################################
  #User Input Section and Setup
  ##################################
  # set paths
  path_dir = '/mnt/bigdata/tetarenk/timing_test/'#make sure to including trailing / !!!
  path_dir_ac='/home/ubuntu/'

  sys.path.append(os.path.join(path_dir_ac, "AstroCompute_Scripts/"))
  ##################################

  ##################################
  #Reading in Parameters
  ##################################
  # Label for casa output directories and files.
  fits_file = path_dir+'data_products/images_Test_Date_6GHz_0hours0min10sec/Test_6GHz_Date_whole_dataset.fits'
  out_file0 = path_dir+'data_products/images_Test_Date_6GHz_0hours0min10sec/Test_6GHz_Date_whole_dataset_aegean.txt'
  tab_file = path_dir+'data_products/images_Test_Date_6GHz_0hours0min10sec/Test_6GHz_Date_whole_dataset_objdet_comp.tab'
  catalog_input_name = path_dir+'data_products/images_Test_Date_6GHz_0hours0min10sec/Test_6GHz_Date_whole_dataset_objdet.tab'
  # aegean parameters
  seed = 10
  flood = 4
  tele = 'VLA'
  cellSize_string = '0.3arcsec'
  ##################################

  #run aegean
  src_l, ra_l, dec_l, maj_l, min_l, pos_l=objdet(tele,lat,out_file0,fits_file,seed,flood,tab_file,catalog_input_name,cellSize_string)

