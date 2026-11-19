# BOUTExtensions

[Root](../../README.md) | [common](../README.md)

Own extensions to `BOUT++`, in C++.

* [include](include) - The `.hxx` files.
* [src](src) - The `.cxx` files.
  Each header ends by including its `.cxx` file, so the models include the
  headers and no library is built from `src`.
* [Doxyfile](Doxyfile) - Doxygen configuration.

Both folders contain the files

* `helpers` - `polAvg` for poloidal averaging and `volumeIntegral` for volume
  integrals.
* `ownBCs` - Boundary conditions, including "boundaries" that treat the
  singularity at the cylinder axis.
* `ownFilters` - Spectral filters.
* `ownLaplacianInversions` - The `NaulinSolver` Laplace inversion.
* `ownMonitors` - Monitors such as the energy and the particle number.
* `ownOperators` - Differencing operators such as `div_f_GradPerp_g` and
  `Grad_perp`. The model equations depend on these; see
  [../../AGENTS.md](../../AGENTS.md).
* `parameters` - A class that converts physical units to normalized
  quantities and calculates the collisionalities and viscosities.

`include` also holds `signedJacobian.hxx`, which has no `.cxx` file.

The operators and boundary conditions are verified by the cases in
[../../MES](../../MES/README.md).
