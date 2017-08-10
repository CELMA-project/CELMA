# celma

The folder containing the `celma` code.
**NOTE**: If the run exits prematurely, see
[../common/CELMAPy/scripts/prematureExitFixes](../common/CELMAPy/scripts/prematureExitFixes)
for common fixes.

## The model
* [CSDXMagFieldScanAr](CSDXMagFieldScanAr) - Contains the `BOUT.inp` file used
  for the magnetic field scan.
* [CSDXNeutralScanAr](CSDXNeutralScanAr) - Contains the `BOUT.inp` file used
  for the neutral density scan.
* [CSDXNyScan](CSDXNyScan) - Contains the `BOUT.inp` file used for the `ny`
  scan.
* [celma.cxx](celma.cxx) - The source file for the `celma` code.
* [celma.hxx](celma.hxx) - The include file for the `celma` code.
* [makefile](makefile) - The `makefile` for the `celma` code.
  **NOTE:** The makefile assumes that `BOUT-dev` is located in `$HOME` through
  the `BOUT-TOP` flag.

## Run scripts
These scripts are made for doing the simulations on a super computer:
* [PBSScanMarconi-CSDXMagFieldScanAr.py](PBSScanMarconi-CSDXMagFieldScanAr.py) - The magnetic field scan.
* [PBSScanMarconi-CSDXNeutralScanAr.py](PBSScanMarconi-CSDXNeutralScanAr.py) - The neutral density scan.
* [PBSScanMarconi-CSDXNyScan.py](PBSScanMarconi-CSDXNyScan.py) - The `ny` scan.

## Post processing scripts
These scripts are made for doing the post processing on a super computer:
* [PBSPlotMarconi-CSDXMagFieldScanAr.py](PBSPlotMarconi-CSDXMagFieldScanAr.py) - To be run when
  [PBSScanMarconi-CSDXMagFieldScanAr.py](PBSScanMarconi-CSDXMagFieldScanAr.py)
  is done.
* [PBSPlotMarconi-CSDXNeutralScanAr.py](PBSPlotMarconi-CSDXMagFieldScanAr.py) - To be run when
  [PBSScanMarconi-CSDXNeutralScanAr.py](PBSScanMarconi-CSDXNeutralScanAr.py)
  is done.
* [PBSPlotMarconi-CSDXNyScan.py](PBSPlotMarconi-CSDXMagFieldScanAr.py) - To be run when
  [PBSScanMarconi-CSDXNyScan.py](PBSScanMarconi-CSDXNyScan.py)
  is done.

## Miscellaneous
* [pickleTweaks](pickleTweaks) - Contains python-scripts which alters
  standard-plots obtained from `PBSPlot<cluster>-<scan>.py`. The scripts were
  used to create the plots given in the
  [thesis](https://github.com/CELMA-project/dissertation/releases/latest).
  **NOTE:** These are to be run *after* the simulations are done.
* [testPostProcessing.py](testPostProcessing.py) - To be run after
  [testScan.py](testScan.py). Checks that the system has been sucessfully built
  for the included post-processing.
* [testScan.py](testScan.py) - Performs a quick scan to check if the system has
  been sucessfully built.
