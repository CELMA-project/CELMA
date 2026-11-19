# celma

[Root](../README.md)

The CELMA model: a `BOUT++` model of plasma turbulence in a cylinder.

## The model

* [celma.cxx](celma.cxx) and [celma.hxx](celma.hxx) - The source of the model.
  The equations are here and in
  [../common/BOUTExtensions](../common/BOUTExtensions/README.md); edit with
  care, see [AGENTS.md](../AGENTS.md).
* [CMakeLists.txt](CMakeLists.txt) - Builds the `celma` executable next to the
  source, as part of the [root build](../README.md#build-and-test).
* [CSDXMagFieldScanAr](CSDXMagFieldScanAr) - The `BOUT.inp` for the magnetic
  field scan, and `editDmpFoldersDict.py`.
* [CSDXNeutralScanAr](CSDXNeutralScanAr) - The `BOUT.inp` for the neutral
  density scan.
* [CSDXNyScan](CSDXNyScan) - The `BOUT.inp` for the `ny` scan. The regression
  test runs the model with this input.

Run a case from this directory with `./celma -d CSDXNyScan`.

The variant with the Boussinesq approximation is in
[../celmaWithBoussinesqApprox](../celmaWithBoussinesqApprox/README.md).

## Python scripts

These were written for the earlier makefile build and use `bout_runners` and
the `CELMAPy` package, which is in a separate repository (CELMA-py) and not
installed by the CI image.
They are not run by the continuous integration.

* [PBSScanMarconi-CSDXMagFieldScanAr.py](PBSScanMarconi-CSDXMagFieldScanAr.py),
  [PBSScanMarconi-CSDXNeutralScanAr.py](PBSScanMarconi-CSDXNeutralScanAr.py)
  and [PBSScanMarconi-CSDXNyScan.py](PBSScanMarconi-CSDXNyScan.py) - The
  scans, for submission to a PBS queue on a super computer.
* [PBSPlotMarconi-CSDXMagFieldScanAr.py](PBSPlotMarconi-CSDXMagFieldScanAr.py),
  [PBSPlotMarconi-CSDXNeutralScanAr.py](PBSPlotMarconi-CSDXNeutralScanAr.py)
  and [PBSPlotMarconi-CSDXNyScan.py](PBSPlotMarconi-CSDXNyScan.py) - The post
  processing of the matching scan, to be run when it is done.
* [pickleTweaks](pickleTweaks) - Scripts that alter the standard plots from
  the post processing, used for the plots in the
  [thesis](https://github.com/CELMA-project/dissertation/releases/latest).
  Run them after the simulations are done.
* [testScan.py](testScan.py) - A quick scan to check that the system is built.
* [testPostProcessing.py](testPostProcessing.py) - Run after `testScan.py`;
  checks the post processing.
* [brokenExitExample.py](brokenExitExample.py) - Example of how to repair
  runs whose `*.dmp.*` files were corrupted by a premature exit while
  writing. Do not run it unless there is a problem.
