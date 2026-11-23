# MES

[Root](../README.md)

Cases that check that the operators and boundary conditions in
[common/BOUTExtensions](../common/BOUTExtensions/README.md) converge at the
expected order.
MES is the method of manufactured solutions: pick an analytic solution,
compute the source it needs, run the operator on a sequence of finer grids
and check that the error falls at the expected order.

Each case directory holds:

* the C++ source and a `CMakeLists.txt` (listed in the
  [root CMakeLists.txt](../CMakeLists.txt));
* `calculations/`: the notebook (`exactSolutions.ipynb`) that derives the
  analytic solution and source, and `functions.txt` where present;
* one subdirectory per analytic solution, holding its `BOUT.inp`;
* `driver*.py` per solution and `conclusion.txt` with the recorded result.
  The drivers use `bout_runners` and `CELMAPy`, which are not in this
  repository, see [AGENTS.md](../AGENTS.md).

[../tests/regression/mes.py](../tests/regression/mes.py) runs every
`BOUT.inp` below at three small resolutions, compares the largest absolute
error to a stored reference and checks the observed order of convergence.
[../tests/regression/known-orders.json](../tests/regression/known-orders.json)
lists the cases that fall short of the order of their scheme, each with its
cause.
The convergence results in the `conclusion.txt` files are from the original
runs and are not rechecked by the continuous integration.

## arakawaBracket

* [arakawaBracket](arakawaBracket)
  * [ArakawaOfDDXf2](arakawaBracket/ArakawaOfDDXf2/conclusion.txt):
      the `bracket(a, b, bm)` method on `{(d_rho phi)^2, n}`. Converges.
      Solution: `mixModeAndGaussian`.

## boundaries

* [boundaries](boundaries)
  * [1-yExtrapolation](boundaries/1-yExtrapolation/conclusion.txt):
      extrapolation to the y ghost points. Fourth order, reaches machine
      precision. Solution: `gaussianWSinAndParabola`.
  * [2-uEParSheath](boundaries/2-uEParSheath/conclusion.txt):
      the parallel electron velocity sheath condition. Converges.
      Solutions: `gaussianWSinAndParabola`, `noYDepInProf`.
  * [3-cauchyBC](boundaries/3-cauchyBC/conclusion.txt):
      the Cauchy boundary condition. Only first order in every trial.
      Solutions: `gaussianWSinAndParabola`, `noXDepInFunc`.
  * [4-jParSheath](boundaries/4-jParSheath/conclusion.txt):
      the parallel current sheath condition. Converges.
      Solution: `gaussianWSinAndParabola`.

## divOfScalarTimesVector

* [divOfScalarTimesVector](divOfScalarTimesVector)
  * [1a-divPerp](divOfScalarTimesVector/1a-divPerp/conclusion.txt):
      `div(a*b)` with a `Vector3D` `b`, perpendicular part. First order; the
      conclusion file points to a `subConvergence.ipynb` that is not in the
      repository. Solution: `twoGaussians`.
  * [1b-JTimesDivPerp](divOfScalarTimesVector/1b-JTimesDivPerp/conclusion.txt):
      the same multiplied by the Jacobian. Second order.
      Solution: `twoGaussians`.
  * [2a-divSource](divOfScalarTimesVector/2a-divSource/conclusion.txt):
      the `divSource` term. Converges for two solutions and
      diverges for `twoNormalGaussians`, cause unknown.
      Solutions: `twoGaussians`, `twoGaussiansSwapped`,
      `twoNormalGaussians`.
  * [2b-JTimesDivSource](divOfScalarTimesVector/2b-JTimesDivSource/conclusion.txt):
      the same multiplied by the Jacobian. Same pattern as 2a.
      Solutions: `twoGaussians`, `twoGaussiansSwapped`,
      `twoNormalGaussians`.

## integrals

* [integrals](integrals)
  * [volumeIntegral](integrals/volumeIntegral/conclusion.txt):
      the volume integral. Converges in x and y; the periodic z direction
      with sines gives no convergence at small errors.
      Solutions: `sumSines`, `sumSinesHat`, `xHat`, `yHat`, `zHat`,
      `zSines`.

## laplaceInversion

* [laplaceInversion](laplaceInversion)
  * [Naulinsolver0Bndry](laplaceInversion/Naulinsolver0Bndry/conclusion.txt):
      inversion of `div(a*Grad_perp(b))-f` with the Naulin solver and zero
      potential on the boundary. Second order.
      Solutions: `twoGaussians`, `twoGaussiansSwapped`,
      `twoNormalGaussians`.

## polAvg

* [polAvg](polAvg/conclusion.txt): the poloidal average. Second order, from
  the rectangular integration. Solution: `zHat`. The conclusion file says
  "To be run".

## singleOperators

Single operators in cylinder geometry, checked for convergence despite the
singularity at the axis.
Two notebooks at this level explain the failing cases:
[subConvergenceDDXOverJ.ipynb](singleOperators/subConvergenceDDXOverJ.ipynb)
and [properZFailConvergence.ipynb](singleOperators/properZFailConvergence.ipynb).

* [singleOperators](singleOperators)
  * [1-DDXCylinder](singleOperators/1-DDXCylinder/conclusion.txt):
      `DDX(f)`. Converges. Solutions: `properZ`, `gaussian`,
      `gaussianWSinAndParabola`.
  * [1.1-DDXCylinderWideStencilCenter](singleOperators/1.1-DDXCylinderWideStencilCenter/conclusion.txt):
      `DDX(f)` with a sixth order stencil at the two innermost points.
      Converges. Solutions: `gaussian`, `gaussianWSinAndParabola`.
  * [2-D2DX2Cylinder](singleOperators/2-D2DX2Cylinder/conclusion.txt):
      `DDX(DDX(f))`. Converges, except `properZ`, which diverges for an
      unknown reason. Solutions: `properZ`, `gaussian`,
      `gaussianWSinAndParabola`.
  * [3-DDXOverJCylinder](singleOperators/3-DDXOverJCylinder/conclusion.txt):
      `(1/J)*DDX(f)`. First order at best, broken convergence for three
      solutions, see the two notebooks above. Solutions: `gaussian`,
      `gaussianWSinAndParabola`, `properZ`, `properZModified`, `simpleX`.
  * [4-DDZCylinder](singleOperators/4-DDZCylinder/conclusion.txt):
      own `DDZ`. Errors at machine precision.
      Solution: `gaussianWSinAndParabola`.
  * [5-D2DZ2Cylinder](singleOperators/5-D2DZ2Cylinder/conclusion.txt):
      own `D2DZ2`. Errors at machine precision.
      Solution: `gaussianWSinAndParabola`.
