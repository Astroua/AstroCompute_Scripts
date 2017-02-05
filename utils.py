
import json
import re
import astroML.time_series
from scipy.stats import norm
import matplotlib.pyplot as pp
import numpy as np
from scipy.stats import chi2
import casac
from tasks import *
from taskinit import *


def convert_param_format(filename, to="json"):
    '''
    Convert the text format into a dictionary, or JSON.
    '''

    with open(filename, "r") as f:
        params = f.read()

    params = params.split("\n")
    param_dict = {}

    for param in params:
        # Skip comment lines
        if len(param) == 0:
            continue

        if param[0] == "#":
            continue

        # Strip whitespace
        param = re.sub("\s", "", param)

        splits = param.split("=")

        key = splits[0]
        value = splits[1]

        if len(value.split("#")) > 1:
            value = value.split("#")[0]

        param_dict[key] = value.strip("'")

    if to is "json":
        return json.dumps(param_dict)
    else:
        return param_dict


def load_json(filename):
    '''
    Load in a JSON formatted text file.
    '''

    with open(filename, "r") as f:
        contents = json.load(f)

    for key in contents:
        if isinstance(contents[key], list):
            continue
        try:
            contents[key] = int(contents[key])
        except ValueError:
            contents[key] = str(contents[key])

    return contents


def is_power2(num):
    ''' Check if input imsize is a power of 2^n, in order to
    optimize cleaning speed

    num: image size in pixels

    return: True/False
    '''
    return(num != 0 and ((num & (num-1)) == 0))

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

def initial_clean(visibility,outputPath,label,imageSize,cellSize,spw_choice,taylorTerms,numberIters,thre):
    clean(vis=visibility,
          imagename=os.path.join(outputPath, label+'whole_dataset'),
          field='', mode='mfs', imsize=imageSize, cell=cellSize,
          weighting='natural', spw=spw_choice, nterms=taylorTerms,
          niter=numberIters, gain=0.1,
          threshold=thre, interactive=False)
    exportfits(imagename=os.path.join(outputPath, label+'whole_dataset.image'),
               fitsimage=os.path.join(outputPath, label+'whole_dataset.fits'),
               history=False)

def errf(ampl,y,er):
   ''' Residual function for Chi2 calculation
   
   ampl: weighted mean
   y: flux array
   er: flux error array
   
   return: residual for chi^2 calculation in chi2_calc
   '''
   fitf = ampl
   return (y-fitf)/er

def chi2_calc(flux,fluxerr):
   ''' Chi2 with constant flux model
   
   flux: flux array
   fluxerr: flux error array
   
   return: chi^2 with constant flux (at weighted mean) model
   '''
   we_fix=[]
   for item in fluxerr:
      w_fix=1/((item)**2)
      we_fix.append(w_fix)
   wei_fix=np.array(we_fix)
   dof_fix=len(flux)-1
   wm_fix=np.average(flux,weights=wei_fix)
   un_fix=1/(np.array(we_fix).sum())
   residual_fix=errf(wm_fix,flux,fluxerr)
   chisquared_fix=residual_fix**2  
   chi_tot_fix=((residual_fix**2).sum())
   null_hyp_fix=chi2.sf(chi_tot_fix,(np.array(flux).shape[0])-1)
   return(chi_tot_fix,dof_fix,wm_fix,un_fix,null_hyp_fix)
	
def lomb_scargle(time,flux,fluxerr,interval,label):
   '''Generalized LS periodogram (Note: Power is normalized between 0 and 1)
   
   time: MJD array
   flux: flux array
   fluxerr: flux error array
   interval: time bin size in seconds
   label: name for savefig
   
   return: plot is saved to a file
   '''

   secondsElapsed=[]
   for i in range(0,len(time)):
      secondsElapsed.append((time[i]-time[0])*24*60*60+interval/2.0)
   omega=np.logspace(np.log10(2.*np.pi/(secondsElapsed[-1])),np.log10(2*np.pi/(2.*interval)),10000)
   samp=1./interval
   lsg,sig=astroML.time_series.lomb_scargle(secondsElapsed,flux,fluxerr,omega,generalized=True,significance=[0.05,0.01])
   fig=pp.figure()
   ax1=fig.add_subplot(111)
   ax1.plot(omega/(2*np.pi),(lsg))
   pp.axhline(y=sig[0],linewidth=4,ls='--',color='m')
   pp.axhline(y=sig[1],linewidth=4,ls='--',color='c')
   pp.xlim(min(omega)/(2.*np.pi),max(omega)/(2.*np.pi))
   #pp.xscale("log")
   #pp.yscale("log")
   pp.xlabel('Frequency, $\\nu$ (Hz)',size=16)
   pp.ylabel('Lomb-Scargle Power',size=16)
   pp.savefig(label)
   return(sig[0],sig[1])

def var_analysis(flux,fluxerr):
   '''Run all variability analysis

   flux: flux array
   fluxerr: flux error array

   return: total chi^2, degrees of freedom, null hypothesis probability, weighted mean,
   weighted mean error, excess variance, excess variance error, fractional rms, fractional rms error
   '''

#chi2 and weighted mean
   chi_tot,dof,wm,wmerr,null=chi2_calc(flux,fluxerr)
#excess variance and fractional rms
   var_data=np.var(flux,ddof=1)
   rms_mean=np.sum(fluxerr**2)/len(fluxerr)
   ex_var=(var_data)-rms_mean
   if ex_var < 0.0:
      ex_var='n/a'
      print 'Variance of data much less then measurment errors'
      frac_rms='n/a'
      frac_rms_err='n/a'
      ex_var_err='n/a'
   else:
      frac_rms=np.sqrt(ex_var/wm**2)
      ex_var_err=np.sqrt((np.sqrt(2/len(flux))*rms_mean/wm**2)**2+(np.sqrt(rms_mean/len(flux))*2.*frac_rms/wm)**2)
      frac_rms_err=(1./2.*frac_rms)*ex_var_err
   return(chi_tot,dof,null,wm,wmerr,ex_var,ex_var_err,frac_rms,frac_rms_err)
