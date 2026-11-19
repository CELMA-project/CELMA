# celmaWithBoussinesqApprox

[Root](../README.md)

The CELMA model with the Boussinesq approximation.
The model without it is in [../celma](../celma/README.md).

## The model

* [celmaWBA.cxx](celmaWBA.cxx) and [celmaWBA.hxx](celmaWBA.hxx) - The source
  of the model. The equations are here and in
  [../common/BOUTExtensions](../common/BOUTExtensions/README.md); edit with
  care, see [AGENTS.md](../AGENTS.md).
* [CMakeLists.txt](CMakeLists.txt) - Builds the `celmaWBA` executable next to
  the source, as part of the [root build](../README.md#build-and-test).
* [BousCSDXMagFieldScanAr](BousCSDXMagFieldScanAr) - The `BOUT.inp` for the
  magnetic field scan. The regression test runs the model with this input.
* [BousCSDXNeutralScanAr](BousCSDXNeutralScanAr) - The `BOUT.inp` for the
  neutral density scan.

Run a case from this directory with `./celmaWBA -d BousCSDXMagFieldScanAr`.

## Python scripts

These were written for the earlier makefile build and use `bout_runners` and
the `CELMAPy` package, which is in a separate repository (CELMA-py) and not
installed by the CI image.
They are not run by the continuous integration.

* [PBSScanMarconi-BousCSDXMagFieldScanAr.py](PBSScanMarconi-BousCSDXMagFieldScanAr.py)
  and [PBSScanMarconi-BousCSDXNeutralScanAr.py](PBSScanMarconi-BousCSDXNeutralScanAr.py) -
  The magnetic field scan and the neutral density scan, for submission to a
  PBS queue.
* [PBSPlotMarconi-BousCSDXMagFieldScanAr.py](PBSPlotMarconi-BousCSDXMagFieldScanAr.py) -
  Post processing of the magnetic field scan, to be run when it is done.
* [testScan.py](testScan.py) - A quick scan to check that the system is built.
* [testPostProcessing.py](testPostProcessing.py) - Run after `testScan.py`;
  checks the post processing.
