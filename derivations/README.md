# derivations

[Root](../README.md)

Derivations of stencils, boundary conditions and other analytical
calculations, mostly in jupyter notebooks.
Nothing here is run by the simulations: the results are copied into the C++ by
hand, so changing a notebook does not change the code.
The notebooks are committed without output.

* [boundaries](boundaries/README.md) - Polynomials for the Cauchy boundary
  condition and for extrapolation to the ghost point.
    * [cauchyBC](boundaries/cauchyBC/README.md)
    * [ghostPointExtrapolation](boundaries/ghostPointExtrapolation/README.md)
* [collisionalities](collisionalities/README.md) - Scripts for the
  electron-neutral and ion-neutral collision frequencies in argon and for
  plasma parameters.
* [dispersion](dispersion/README.md) - Dispersion relation for slab geometry
  without simplifications.
* [divOfExBOperator](divOfExBOperator/README.md) - The divergence of the
  `E x B` advection of the perpendicular gradient of the potential times the
  density, in Clebsch coordinates.
* [laplacianInversion](laplacianInversion/README.md) - Stencils for
  `div(A*Grad_perp(f))+Bf`, kept for a possible non-Boussinesq Laplace
  inversion.

The cases that check the resulting operators are in
[../MES](../MES/README.md).
