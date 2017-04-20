# BOUTExtensions

Contains own extensions to `BOUT++`.

* [include](include) - Contains the `.hxx` files.
* [src](src) - Contains the `.cxx` files.

Both folders contains the files

* `helpers` - Contains `polAvg` for poloidal averaging, and `volumeIntegral`
  for volume integrals
* `ownBCs` - Contains boundary conditions.
  Specifically:
    * "Boundaries" treating the singularity at the cylinder axis.
* `ownFilters` - Contain different spectral filters.
* `ownLaplacianInversions` - Contains the `NaulinSolver` laplacian inversion .
* `ownMonitors` - Contains monitors like the energy and particle number.
* `ownOperators` - Contains differencing operators like `div_f_GradPerp_g` and
  `Grad_perp`.
* `parameters` - Contains a class which handles the conversion from physical
  units to normalized quantites. This class also calculates the collisionalites
  and viscosities.
