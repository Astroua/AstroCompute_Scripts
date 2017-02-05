# AstroCompute CASA Scripts
A collection of python scripts to create high time resolution light curves from calibrated interferometric data sets. These scripts run in the Common Astronomy Software Application ([CASA](http://casa.nrao.edu)) and have been tested with VLA, SMA and NOEMA data sets.

## Links for importing data into CASA
* VLA: Can be downloaded from archive in CASA format.
* SMA: Can be calibrated in CASA. Follow instructions [here](https://www.cfa.harvard.edu/sma/casa) for reduction details and details on how to convert to a CASA MS data set for use by these scripts
* NOEMA: Must be calibrated beforehand. Follow instructions [here](http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf) to convert to a CASA MS data set for use by these scripts.

##Requires the following python packages
* casa-python executable wrapper (get it [here](https://github.com/radio-astro-tools/casa-python))
   * use this to install these python packages within CASA python
      * jdcal (`casa-pip install jdcal`)
      * astropy (`casa-pip install astropy`)
      * astroML (`casa-pip install astroML`)
* uvmultifit (get it [here](http://nordic-alma.se/support/software-tools))
   * this package also needs gsl libraries to build (get them [here](http://askubuntu.com/questions/490465/install-gnu-scientific-library-gsl-on-ubuntu-14-04-via-terminal))
* analysisUtils (get it [here](https://casaguides.nrao.edu/index.php?title=Analysis_Utilities))
* (*optional*) aegean (get it [here](https://github.com/PaulHancock/Aegean))
   * this package also needs lmfit-0.7.4 (get it [here](http://github.com/lmfit/lmfit-py.git@0.7.4), install with `casa-pip install git+git://github.com/lmfit/lmfit-py.git@0.7.4`)


##Description of scripts
1. Primary script:
   * **casa_timing_script.py**: intended to be run within CASA. This is the script that does all the hard work.
      * All parameters need to be carefully considered and changed for each new data set.
      * If you have a complicated field with other sources it is recommended that you use your own mask file (with clean boxes     around bright sources; mask_option='file') for cleaning, or run object detection, which will create a mask file with a       boxed region around each detected source (mask_option='aegean').
      * Included in casa_timing_script.py is the option to run basic variability analysis:
         * Calculate weighted mean and do a chi^2 with a constant flux model,
         * Calculate excess variance (Vaughn et al. 2003; Bell et al., 2015),
         * Calculate fractional RMS (Vaughn et al. 2003; Bell et al., 2015), and
         * Make power spectrum using generalized lomb-periodogram (Zechmeister and Kurster, 2009)
2. Other scripts:
   * **Aegean_ObjDet.py**: object detection algorithm.
      * Is integrated into casa_timing_script.py, but can be run on its own too.
      * A data file of detected object properties is output.
   * **utils.py**: module of tools used in casa_timing_script.py.
   * **download_data_AWS.py**: downloads data from an AWS bucket.
   * **upload_data_AWS.py**: uploads data to an AWS bucket.
   * **remove_bucket_AWS.py**: removes a bucket from AWS.
   * **sim.py**: creates a simulated CASA MS of time-variable source for testing using CASA's simulation toolbox.

##User Parameters
* All input parameters for casa_timing_script.py are set in the parameter file (param.txt)
* A complete description of each parameter is provided in param_description.txt.

####
Support from the SKA/AWS AstroCompute in the Cloud Program
