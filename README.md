# AstroCompute CASA Scripts
A collection of python scripts to create high time resolution light curves from calibrated interferometric data sets. These scripts run in the Common Astronomy Software Application ([CASA](http://casa.nrao.edu)) and have been tested with VLA, ALMA, SMA, ATCA, and NOEMA data sets.

## Importing data into CASA
* VLA: Can be directly imported into CASA, and downloaded from the archive in CASA format.
* ALMA: Can be directly imported into CASA, and downloaded from the archive in CASA format.
* SMA: Can be imported into CASA with a little extra work. Follow instructions [here](https://www.cfa.harvard.edu/sma/casa) on converting to a CASA MS data set.
* NOEMA: Must be calibrated beforehand, but can be imaged in CASA. Follow instructions [here](http://www.iram.fr/IRAMFR/ARC/documents/filler/casa-gildas.pdf) for converting to a CASA MS data set.
* ATCA: Can be directly imported into CASA. Follow instructions [here](https://casaguides.nrao.edu/index.php/ATCA_Tutorials) on converting to a CASA MS data set.

## Requires
* **CASA** (get it [here](https://casa.nrao.edu/casa_obtaining.shtml), scripts are tested up to version 5.4)
* The following **python packages** need to be installed within CASA; **jdcal, astropy, astroML**<br/>
For importing python packages into CASA ([detailed instructions](http://docs.astropy.org/en/stable/install.html)),<br/>
```
casa --no-logger --log2term -c "from setuptools.command import easy_install; easy_install.main(['--user', 'pip'])"
casa --no-logger --log2term -c "import pip; pip.main(['install', 'astropy', '--user']); pip.main(['install', 'astroML', '--user']); pip.main(['install', 'jdcal', '--user'])"
```
* **analysisUtils** (get it [here](https://casaguides.nrao.edu/index.php?title=Analysis_Utilities))
* (*optional*) For UV plane fitting, **uvmultifit** (get it [here](http://nordic-alma.se/support/software-tools))
* (*optional*) To use object detection, you need the **aegean** package (see [here](https://github.com/PaulHancock/Aegean)) <br/>
```
casa --no-logger --log2term -c "import pip; pip.main(['install', 'git+https://github.com/PaulHancock/Aegean.git', '--user'])"
```

## To make your own high time-resolution lightcurves you will need:
1. **casa_timing_script.py**: intended to be run within CASA. This is the script that does all the hard work.
2. **param.txt**: All input parameters are set in this parameter file. A complete description of each parameter is provided in `param_description.txt`.
2. **utils.py**: module of tools.
3. (*optional*), **Aegean_ObjDet.py**: object detection algorithm.
   * Integrated into the timing script, but can also be run on its own, with a fits image as input. Output is a labelled image and region files (DS9 and CASA format) of detected sources.
## How to use the scripts:
```
casa -c casa_timing_script.py [path_to_param_file] [path_dir] [path_to_repo]
```
   * All MSs need to be in `path_dir/data`, all output goes into `path_dir/data_products` (directory created by script).
   * Make sure to include trailing `/` in all directory names.
   * Please create a split MS, with only the target present.
   * All parameters need to be carefully considered and changed for each new data set.
   * If you have a complicated field with other sources it is recommended that you use your own mask file (with clean boxes around bright sources; `mask_option='file'`) for cleaning, or run object detection, which will create a mask file with a boxed region around each detected source (`mask_option='aegean'`).
   * If you are using an outlier field file, please follow the example template in `example_outlier_file.txt`.
   * If you are doing uv-fitting, please follow the example initial parameters template in `uv_init_example.txt`.
   * The option to run basic variability analysis is included (`var_anal='T'`):
      * Calculate weighted mean and do a chi^2 with a constant flux model
      * Calculate excess variance and fractional RMS ([Vaughan et al. 2003](http://adsabs.harvard.edu/abs/2003MNRAS.345.1271V); [Bell et al., 2015](http://adsabs.harvard.edu/abs/2015MNRAS.450.4221B))
      * Make a power spectrum using the generalized lomb-periodogram algorithm ([Zechmeister and Kurster, 2009](http://adsabs.harvard.edu/abs/2009A%26A...496..577Z))

####
Support from the SKA/AWS AstroCompute in the Cloud Program
