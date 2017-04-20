# MES

This folder contains checks to see that the implementation is convergent.

* [arakawaBracket](arakawaBracket) - Verifies the `bracket(a, b, bm)` method.
* [boundaries](boundaries) - Verifies the boundary conditions given in
  `ownBCs`.
* [divOfScalarTimesVector](divOfScalarTimesVector) - Verifies terms on the form
  `div(a*b)`, where `b` is a `Vector3D` object.
* [integrals](integrals) - Verifies the integrations shcemes.
* [laplaceInversion](laplaceInversion) - Verifies the inversion of
  `div(a*Grad_perp(b))-f`.
* [polAvg](polAvg) - Verifies the poloidal average.
* [singleOperators](singleOperators) - Verifies single operators like `DDX(a)`,
  to see that convergence is given despite the singularity at the axis.
