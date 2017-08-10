# celmaWithBoussinesqApprox

The folder containing the `celma` code, where the Boussinesq
approximation has been used.
**NOTE**: If the run exits prematurely, see
[../common/CELMAPy/scripts/prematureExitFixes](../common/CELMAPy/scripts/prematureExitFixes)
for common fixes.

## The model
* [BousCSDXMagFieldScanAr](BousCSDXMagFieldScanAr) - Contains the
  `BOUT.inp` file used for the magnetic field scan.
* [BousCSDXNeutralScanAr](BousCSDXNeutralScanAr) - Contains the
  `BOUT.inp` file used for the neutral density scan.
* [celmaWBA.cxx](celmaWBA.cxx) - The source file for the "`celma` with
  Boussinesq Approximation" code.
* [celmaWBA.hxx](celmaWBA.hxx) - The include file for the "`celma` with
  Boussinesq Approximation" code.
* [makefile](makefile) - The `makefile` the "`celma` with
  Boussinesq Approximation" code.
  **NOTE:** The makefile assumes that `BOUT-dev` is located in `$HOME` through
  the `BOUT-TOP` flag.

## Run scripts
These scripts are made for doing the simulations on a super computer:
* [PBSPlotMarconi-BousCSDXMagFieldScanAr.py](PBSScanMarconi-CSDXMagFieldScanAr.py) -
  The magnetic field scan for the Boussinesq approximated `celma` model.
* [PBSScanMarconi-BousCSDXMagFieldScanAr.py](PBSScanMarconi-BousCSDXMagFieldScanAr.py) -
  The neutral density scan for the Boussinesq approximated `celma` model.

## Post processing scripts
These scripts are made for doing the post processing on a super computer:
* [PBSPlotMarconi-BousCSDXMagFieldScanAr.py](PBSPlotMarconi-BousCSDXMagFieldScanAr.py) - To be run when
  [PBSScanMarconi-BousCSDXMagFieldScanAr.py](PBSScanMarconi-BousCSDXMagFieldScanAr.py)
  is done.

## Miscellaneous
* [pickleTweaks](pickleTweaks) - Contains python-scripts which alters
  standard-plots obtained from `PBSPlot<cluster>-<scan>.py`. The scripts were
  used to create the plots given in the
  [thesis](https://github.com/CELMA-project/dissertation/releases/latest).
  **NOTE:** These are to be run *after* the simulations are done.
* [brokenExitExample.py](brokenExitExample.py) - Example on how to repair
  simulations where the `*.dmp.*`-files have been corrupted due to for example
  premature exit during the flush stage of writing.
* [testPostProcessing.py](testPostProcessing.py) - To be run after
  [testScan.py](testScan.py). Checks that the system has been sucessfully built
  for the included post-processing.
* [testScan.py](testScan.py) - Performs a quick scan to check if the system has
  been sucessfully built.
