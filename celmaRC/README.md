# celmaCurMom

Folder with the different versions of the code when evolving the parallel current and
total parallel density momentum.

Each run of each group has its own driver. `bout_runners` is used for running
and post-processing.

* 1-CELMA-RC - As `celmaCurMom/7-CELMACurMomWParamsConstViscClean`, but without
  the `0.51` bug
* 2-CELMA-RCParHyperVsic - As 1-CELMA-RC, but with hyper viscosities
* 3-CELMA-RCUpwind - As 1-CELMA-RC, but with upwinding on `n*DDY(lnN-phi)`
  term. The problem is always failing in the `Naulin Solver`
* B1-CELMA-RC - As `celmaCurMom/B7-CELMACurMomWParamsConstViscClean`, but without
  the `0.51` bug
* common - python post processing and own implementations to BOUT++
