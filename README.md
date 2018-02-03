# AstroCompute CASA Scripts
A collection of python scripts to create high time resolution light curves from calibrated interferometric data sets. These scripts run in the Common Astronomy Software Application ([CASA](http://casa.nrao.edu)) and have been tested with VLA, SMA and NOEMA data sets.

## Links for importing data into CASA
* VLA: Can be downloaded from archive in CASA format.
* SMA: Can be calibrated in CASA. Follow instructions [here](https://www.cfa.harvard.edu/sma/casa) for reduction details and details on how to convert to a CASA MS data set for use by these scripts
* NOEMA: Must be calibrated beforehand. Follow instructions [here](http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf) to convert to a CASA MS data set for use by these scripts.

## Requires
* CASA (get it [here](https://svn.cv.nrao.edu/casa/linux_distro/release/el6), scripts tested with version 5.1.2)
* The following python packages need to be installed within CASA; `jdcal`, `astropy`, `astroML` <br/>
If using below CASA v5,
   * Use casa-python executable wrapper (get it [here](https://github.com/radio-astro-tools/casa-python)) <br/>
   e.g., `casa-pip install jdcal`
If using CASA v5 and above,
   * Use instructions [here](http://docs.astropy.org/en/stable/install.html) <br/>
   e.g., ```casa --no-logger --log2term -c "from setuptools.command import easy_install; easy_install.main(['--user', 'pip'])"```
   ```casa --no-logger --log2term -c "import pip; pip.main(['install', 'astropy', '--user']); pip.main(['install', 'astroML', '--user']); pip.main(['install', 'jdcal', '--user'])"```
* (*optional*)
   * For UV plane fitting, `uvmultifit` (get it [here](http://nordic-alma.se/support/software-tools))
   * To use object detection,
      * `analysisUtils` (get it [here](https://casaguides.nrao.edu/index.php?title=Analysis_Utilities))
      * `aegean` (see [here](https://github.com/PaulHancock/Aegean)) <br/>
      To install, ```casa --no-logger --log2term -c "import pip; pip.main(['install', 'git+https://github.com/PaulHancock/Aegean.git', '--user'])"```

## Description of scripts
1. Primary script:
   * **casa_timing_script.py**: intended to be run within CASA. This is the script that does all the hard work.
      * All parameters need to be carefully considered and changed for each new data set.
      * If you have a complicated field with other sources it is recommended that you use your own mask file (with clean boxes     around bright sources; mask_option='file') for cleaning, or run object detection, which will create a mask file with a       boxed region around each detected source (mask_option='aegean').
      * If you are using an outlier field file, please follow the example template in example_outlier_file.txt.
      * Included in casa_timing_script.py is the option to run basic variability analysis:
         * Calculate weighted mean and do a chi^2 with a constant flux model,
         * Calculate excess variance,
         * Calculate fractional RMS ([Vaughan et al. 2003](http://adsabs.harvard.edu/abs/2003MNRAS.345.1271V); [Bell et al., 2015](http://adsabs.harvard.edu/abs/2015MNRAS.450.4221B)), and
         * Make power spectrum using generalized lomb-periodogram ([Zechmeister and Kurster, 2009](http://adsabs.harvard.edu/abs/2009A%26A...496..577Z))
2. Other scripts:
   * **Aegean_ObjDet.py**: object detection algorithm.
      * Integrated into casa_timing_script.py.
      * Script can also be run on its own, with a fits image as input. Output is a labelled image and region files (DS9 and CASA format) of detected sources.
   * **utils.py**: module of tools used in casa_timing_script.py.
   * **download_data_AWS.py**: downloads data from an AWS bucket.
   * **upload_data_AWS.py**: uploads data to an AWS bucket.
   * **remove_bucket_AWS.py**: removes a bucket from AWS.

## User Parameters
* All input parameters for casa_timing_script.py are set in the parameter file (param.txt)
* A complete description of each parameter is provided in param_description.txt.

####
Support from the SKA/AWS AstroCompute in the Cloud Program
