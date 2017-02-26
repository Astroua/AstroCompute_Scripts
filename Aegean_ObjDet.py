##################################
#Object Detection Script
##################################
'''Uses AEGEAN algorithm (https://github.com/PaulHancock/Aegean) for object detection in a radio image.
INPUT: Cleaned FITS image of whole data set--> [fits_file]
OUTPUT: (1) Data file of the parameters of objects found in image--> [tab_file]
        (2) CASA region file of detected sources--> [path_dir]_casa_region.txt
        (3) DS9 region file of detected sources--> [path_dir]_ds9_region.reg
        (4) Image of detected sources--> [path_dir]_detected_sources.pdf
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
from astropy import units as u
from astropy.io import fits
import pyfits
import pylab as pl
import math as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors

def objdet(tele,lat,out_file0,fits_file,seed,flood,tab_file,catalog_input_name,cellSize_string):
  sources = []

  # get latitude for telescope
  
  # adding SMA and NOEMA to list
  if tele == 'SMA':
    lat = 19.8243
  elif tele == 'NOEMA':
    lat = 44.6339
  else:
    lat = aegean.scope2lat(tele)
  

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

def casa_reg_file(src_l,ra_l,dec_l,maj_l,min_l,pos_l,filename,wmap1):
  '''Write CASA region file with bounding boxes around each detected source.
  NOTE: Can be used in viewer or as a mask file in CLEAN.
  '''
  mask_file=open(filename,'w')
  mask_file.write('{0}\n'.format('#CRTFv0'))
  mask_file.write('{0}\n'.format('global color=magenta'))
  for i in range(0,len(src_l)):
    ra=ra_l[i]
    dec=dec_l[i]
    ras=ra[0:2]+'h'+ra[3:5]+'m'+ra[6:]+'s'
    decs=dec[0:3]+'d'+dec[4:6]+'m'+dec[7:]+'s'
    coord0=SkyCoord(ras,decs,frame='icrs')
    x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
    y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
    pos=[x1,y1]
    bbox_halfwidth=np.sqrt((min_l[int(i)-1]*np.cos(pos_l[int(i)-1]))**2+(min_l[int(i)-1]*np.sin(pos_l[int(i)-1]))**2)+3
    bbox_halfheight=np.sqrt((maj_l[int(i)-1]*np.cos(pos_l[int(i)-1]+(np.pi/2.)))**2+(maj_l[int(i)-1]*np.sin(pos_l[int(i)-1]+(np.pi/2.)))**2)+3
    Box = str(pos[0]-bbox_halfwidth)+','+ str(pos[1]-bbox_halfheight)+','+str(pos[0]+bbox_halfwidth)+','+ str(pos[1]+bbox_halfheight)
    mask_file.write('{0}\n'.format('box[['+Box.split(',')[0]+'pix,'+Box.split(',')[1]+'pix], ['+Box.split(',')[2]+'pix,'+Box.split(',')[3]+'pix]], label='+'"'+str(src_l[i])+'"'))
  mask_file.close()
def ds9_reg_file(src_l,ra_l,dec_l,maj_l,min_l,pos_l,filename,rad):
  '''Write DS9 region file with bounding ellipses around each detected source
  '''
  reg_file=open(filename,'w')
  reg_file.write('{0}\n'.format('# Region file format: DS9'))
  reg_file.write('{0}\n'.format('global color=magenta'))
  for i in range(0,len(src_l)):
    reg_file.write('{0}\n'.format('icrs;ellipse('+str(ra_l[i])+','+str(dec_l[i])+','+str(rad)+'"'+','+ str(rad)+'"'+','+ '0) # text={'+str(src_l[i])+'} textangle=30'))
  reg_file.close()

if __name__ == "__main__":
  ##################################
  #User Input Section and Setup
  ##################################
  # set paths
  path_dir = '/path/to/dir/'#make sure to including trailing / !!!
  path_dir_ac='/path/to/repo/'

  sys.path.append(os.path.join(path_dir_ac, "AstroCompute_Scripts/"))
  ##################################

  ##################################
  #Reading in Parameters
  ##################################
  # Label for casa output directories and files.
  fits_file = path_dir+'file.fits'
  out_file0 = path_dir+'Source_XXGHz_Date_whole_dataset_aegean.txt'
  tab_file = path_dir+'Source_XXGHz_Date_whole_dataset_objdet_comp.tab'
  catalog_input_name = path_dir+'Source_XXGHz_Date_whole_dataset_objdet.tab'
  # aegean parameters
  seed = 10
  flood = 4
  tele = 'VLA'
  lat=''
  cellSize_string = '0.2arcsec'
  ##################################

  #run aegean
  src_l, ra_l, dec_l, maj_l, min_l, pos_l=objdet(tele,lat,out_file0,fits_file,seed,\
    flood,tab_file,catalog_input_name,cellSize_string)
  #read in fits image file for plotting and get wcs header
  fits_file1=fits_file
  hdulist1 = fits.open(fits_file1)[0]
  data1=hdulist1.data
  wmap1=wcs.WCS(hdulist1.header)
  #write region files of detected sources
  casa_reg_file(src_l,ra_l,dec_l,maj_l,min_l,pos_l,path_dir+'casa_region.txt',wmap1)
  ds9_reg_file(src_l,ra_l,dec_l,maj_l,min_l,pos_l,path_dir+'ds9_region.reg',3)
  #plot map with all detected sources labelled
  fig=plt.figure()
  plt.rc('xtick.major', size=4)
  plt.rc('xtick', color='w', labelsize='large')
  ax1 = fig.add_subplot(111, projection=wmap1.celestial)
  im=plt.imshow(np.nan_to_num(data1[0,0,:,:])*1e6,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=0.5),vmin=0.0,vmax=300.)
  cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
  cbar.set_label('uJy/beam')
  for i in range(0,len(src_l)):
    ra=ra_l[i]
    dec=dec_l[i]
    ras=ra[0:2]+'h'+ra[3:5]+'m'+ra[6:]+'s'
    decs=dec[0:3]+'d'+dec[4:6]+'m'+dec[7:]+'s'
    coord0=SkyCoord(ras,decs,frame='icrs')
    x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
    y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
    e1 = patches.Ellipse((x1,y1), 10, 10,angle=0, linewidth=3, fill=False,color='m')
    ax1.add_patch(e1)
  ax1.text(x1+15,y1+15,str(src_l[i]),color='m',fontsize=15)
  ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
  ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
  ax1.coords['ra'].set_axislabel('Right Ascension')
  ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
  ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
  ax1.set_ylim(0,2048)
  ax1.set_xlim(0,2048)
  plt.savefig(path_dir+'detected_sources.pdf',bbox_inches='tight')
  plt.show()
