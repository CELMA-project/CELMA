# laplacianInversion

Derivations of `b=div(A*Grad_perp(f))+Bf` stencils, which can be used for
non-Boussinesq Laplace inveresion.

At the present time, non such implementation is done in `CELMA`, but the
derivations are left here for possible future use.

* [BackwardsForwards.ipynb](BackwardsForwards.ipynb) - Derives a symmetric
  stencil using a backwards-forwards discretization.
* [ForwardsBackwards.ipynb](ForwardsBackwards.ipynb) - Derives a symmetric
  stencil using a forwards-backwards discretization.
* [ForwardsBackwardsNonSymmetric.ipynb](ForwardsBackwardsNonSymmetric.ipynb) -
  Derives a non-symmetric stencil as the Jacobian is causing a symmetry split.
