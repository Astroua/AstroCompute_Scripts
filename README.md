# AstroCompute CASA Scripts
A collection of python scripts to create high time resolution light curves from calibrated interferometric data sets. These scripts run in the Common Astronomy Software Application (CASA; http://casa.nrao.edu) and have been tested with VLA, SMA and NOEMA data sets.

## Links for importing data into CASA
VLA: Can be downloaded from archive in CASA format.
SMA: Can be calibrated in CASA. Follow instructions at https://www.cfa.harvard.edu/sma/casa for reduction details and details on how to convert to a CASA MS data set for use by these scripts
NOEMA: Must be calibrated beforehand. Follow instructions at http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf to convert to a CASA MS data set for use by these scripts.

###Description of scripts
1. **initial_clean.py**: intended to be run first within CASA to CLEAN whole data set and make a fits image for object detection.
2. **Aegean_ObjDet.py**: intended to be run second outside CASA to detect objects in the whole data set fits image
   * Data file of detected object properties that is read in casa_timing_script.py is output
   * Detected object positions are printed to terminal (target object number needs to be set in parameter file; ind=1)
3. **casa_timing_script.py**: intended to be run within CASA. This script does all the hard work.
   * All input parameters are in the parameter file (param.txt), a complete description of each parameter is provided in param_description.txt.
   * All parameters need to be carefully considered and changed for each new data set.
   * If you have a complicated field with other sources it is recommended that you use your own mask file (with clean boxes     around bright sources; mask_option='file') for cleaning, or run object detection, which will create a mask file with a       boxed region around each detected source (mask_option='aegean').
   * Included in casa_timing_script.py is the option to run basic variability analysis:
      * Calculate weighted mean and do a chi^2 with a constant flux model,
      * Calculate excess variance (Vaughn et al. 2003; Bell et al., 2015),
      * Calculate fractional RMS (Vaughn et al. 2003; Bell et al., 2015), and
      * Make power spectrum using generalized lomb-periodogram (Zechmeister and Kurster, 2009); implemented with
       the astroML package.


####Additional Notes
If you already have a target source position (i.e., target box in pixels, or your own mask image) you can toggle off object detection (skip 1 & 2 and go directly to 3).

####
Support from the SKA/AWS AstroCompute Program
