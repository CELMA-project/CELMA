# laplacianInversion

[Root](../../README.md) | [derivations](../README.md)

Derivations of `b=div(A*Grad_perp(f))+Bf` stencils, which can be used for
non-Boussinesq Laplace inversion.

No such implementation exists in `CELMA` at present; the derivations are left
here for possible future use.

* [BackwardsForwards.ipynb](BackwardsForwards.ipynb) - Derives a symmetric
  stencil using a backwards-forwards discretization.
* [ForwardsBackwards.ipynb](ForwardsBackwards.ipynb) - Derives a symmetric
  stencil using a forwards-backwards discretization.
* [ForwardsBackwardsNonSymmetric.ipynb](ForwardsBackwardsNonSymmetric.ipynb) -
  Derives a non-symmetric stencil, as the Jacobian causes a symmetry split.
