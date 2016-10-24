# celmaCurMom

Folder with the different versions of the code when evolving the parallel current and
total parallel density momentum.

Each run of each group has its own driver. `bout_runners` is used for running
and post-processing.

* 1-CELMACurMom - Original code, copied from `celma/8.3-CELMASplitCleanUp`
* 2.0-CELMACurMomIMEX - As 1, but using IMEX
* 2.1-CELMACurMomIMEXJParEdit - As 2.0, but added stiff jPar terms to the
  diffusive part
* 2.2-CELMACurMomIMEXJParEdit2 - As 2.1, but removed the least stiff jPar terms
  to the diffusive part
* 2.3-CELMACurMomIMEXJParEdit3 - As 2.2, but removed the least stiff jPar terms
  to the diffusive part
* 3.0-CELMACurMomWMonitors - As 1, but added monitors
* 4.0-CELMACurMomWParams - As 2.0, but code input is now done consistently
  through radius, length, n0, B0, Te0, Ti0, and Sn
* 4.1-CELMACurMomWParamsConstVisc - As 4, but viscosity is kept constant
* 5.0-CELMACurMomVisc - As 4.0, but added proper viscosity terms
* 6.0-CELMAHoleOldVisc - As 4.0, but boundary conditions are so that there is a
  "hole" in the sheath where the velocity is 0
* B4.1-CELMACurMomWParamsConstVisc - As 4.1, but using the simple Boussinesq
  approximation
* common - python post processing and own implementations to BOUT++
