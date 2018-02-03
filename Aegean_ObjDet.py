##################################
#Object Detection Script
##################################
'''Uses AEGEAN algorithm (https://github.com/PaulHancock/Aegean) for object detection in a radio image.
INPUT: Cleaned FITS image of whole data set--> [fits_file line 167]
OUTPUT: (1) Data file of the parameters of objects found in image--> [out_file0 line 168]
        (2) CASA region file of detected sources--> casa_region.txt
        (3) DS9 region file of detected sources--> ds9_region.reg
        (4) Image of detected sources--> detected_sources.pdf
NOTE: - Needs lmfit,scipy,astropy pkgs, can now pip install AegeanTools!!
      - Make sure aegean code is in your python path [sys.path.append(path_to_aegean)]

Written by: A. Tetarenko
Last Updated: November 13 2017

TO RUN SCRIPT independently, go to line 155, set variables in User Input Section and Setup &
Reading in Params sections, then run either,
python Aegean_ObjDet.py or casa -c  Aegean_ObjDet.py
'''

# Import modules
import scipy
import lmfit
import astropy
import logging
import logging.config
import sys
import numpy as np
import re
import os
import warnings
import AegeanTools
from AegeanTools.catalogs import save_catalog
from AegeanTools.source_finder import SourceFinder,check_cores
logging.basicConfig(format="%(module)s:%(levelname)s %(message)s")
log = logging.getLogger("Aegean")
# source finding object
sf = SourceFinder(log=log)
from multiprocessing import cpu_count
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

  sf = SourceFinder()

  print 'Running Aegean Object Detection -->'
  detections = sf.find_sources_in_image(fits_file,
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
  sources = sf.sources
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
def run_aegean(tables,cellSize_string):
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
  #cellSize_string=cellSize[0]
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
  # set path for input and output data products
  path_dir = '/path/to/dir/'#make sure to including trailing / !!!
  ##################################

  ##################################
  #Reading in Parameters
  ##################################
  #input fits file name
  fits_file = path_dir+'target_source_image.fits'
  #output detected sources file name
  out_file0 = path_dir+'Source_XXGHz_Date_whole_dataset_aegean.txt'
  ## aegean parameters
  seed = 10
  flood = 4
  tele = 'VLA'
  lat=''#only need if you get error that telescope not in predefined list; in degrees
  cellSize_string = '0.2arcsec'
  ## params for displayed image of detetced sources
  imsize=1024
  imunit='m'#intensity unit in image; 'm'=mJy, 'u'=uJy, 'n'=nJy, ''=Jy
  gam=0.3#powerlaw scaling in image intensity; 1=linear scaling
  vmax=15.#max intensity in image in imunit units
  ##################################

  ##################################
  #Object Detection
  ##################################
  #defines file paths for aegean
  tab_file = out_file0.strip('.txt')+'_objdet_comp.tab'
  catalog_input_name = out_file0.strip('.txt')+'_objdet.tab'
  #run aegean
  src_l, ra_l, dec_l, maj_l, min_l, pos_l=objdet(tele,lat,out_file0,fits_file,seed,\
    flood,tab_file,catalog_input_name,cellSize_string)
  print 'The number of sources detected is: ', len(src_l)
  ##################################

  ##################################
  #Plotting and Region Files
  #of Detected Sources
  ##################################
  if len(src_l)>0:
    print 'Writing region files and plotting labelled map of detected sources...'
    #read in fits image file for plotting and get wcs header
    fits_file1=fits_file
    hdulist1 = fits.open(fits_file1)[0]
    data1=hdulist1.data
    wmap1=wcs.WCS(hdulist1.header)
    #write region files of detected sources
    os.system('rm -rf '+path_dir+'casa_region.txt')
    os.system('rm -rf '+path_dir+'ds9_region.reg')
    casa_reg_file(src_l,ra_l,dec_l,maj_l,min_l,pos_l,path_dir+'casa_region.txt',wmap1)
    ds9_reg_file(src_l,ra_l,dec_l,maj_l,min_l,pos_l,path_dir+'ds9_region.reg',3)
    #set units
    if imunit=='m':
      scaler=1e3
    elif imunit=='u':
      scaler=1e6
    elif imunit=='n':
      scaler=1e9
    elif imunit=='':
      scaler=1.
    #plot map with all detected sources labelled
    fig=plt.figure()
    plt.rc('xtick.major', size=4)
    plt.rc('xtick', color='w', labelsize='large')
    ax1 = fig.add_subplot(111, projection=wmap1.celestial)
    im=plt.imshow(np.nan_to_num(data1[0,0,:,:])*scaler,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=gam),vmin=0.0,vmax=vmax)
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
    ax1.set_ylim(0,imsize)
    ax1.set_xlim(0,imsize)
    plt.savefig(path_dir+'detected_sources.pdf',bbox_inches='tight')
    plt.show()
  else:
    print 'No sources detected, so no region files or plot written.'
  ##################################

  #remove temp files and .last/.log files created by CASA/ipython
  os.system('rm -rf *.last')
  os.system('rm -rf casa*.log')
  print '*********************************************************'
  print 'Script finished. Please inspect resulting data products'
  print '*********************************************************'
  ##################################
