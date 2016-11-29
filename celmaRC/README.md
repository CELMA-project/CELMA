# celmaCurMom

Folder with the different versions of the code when evolving the parallel current and
total parallel density momentum.

Each run of each group has its own driver. `bout_runners` is used for running
and post-processing.

* 1-CELMA-RC - As `celmaCurMom/7-CELMACurMomWParamsConstViscClean`, but without
  the `0.51` bug
* B1-CELMA-RC - As `celmaCurMom/B7-CELMACurMomWParamsConstViscClean`, but without
  the `0.51` bug
* common - python post processing and own implementations to BOUT++
