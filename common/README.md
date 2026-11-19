# common

[Root](../README.md)

Own extensions to `BOUT++`, shared by the models in
[../celma](../celma/README.md) and
[../celmaWithBoussinesqApprox](../celmaWithBoussinesqApprox/README.md) and by
the cases in [../MES](../MES/README.md).

* [BOUTExtensions](BOUTExtensions/README.md) - The C++ extensions: operators,
  boundary conditions, filters, monitors, Laplace inversion and parameters.

The `python` modules for pre- and post-processing (`CELMAPy`) and the
`standardPlots` helpers used to live here.
They moved to a separate repository (CELMA-py), and the Python scripts in
this repository that import them do not run without it.
