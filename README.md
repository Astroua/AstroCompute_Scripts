# AstroCompute_Scripts
CASA timing scripts-->All input parameters in parameter file (param.txt)

1. Run initial_clean.py first within CASA to CLEAN whole data set and make a fits image for object detection
2. Run Aegean_ObjDet.py outside CASA to detect objects in fits image- output is a) data file of detected object properties that is read in timing script and b) object positions are printed to terminal (target object number needs to be set in parameter file; i.e., ind=1)
Note: If you already have a target source position (i.e., target box in pixels) you can toggle off object detection (skip 2 and go directly to 3).
3. Run casa_timing_script.py within CASA
