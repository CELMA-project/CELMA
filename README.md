# CELMA-dev

Developement branch of the CELMA code, copied from revision 1130 of
https://svn.fysik.dtu.dk/mmag/Code-developement/BOUT/

Current implementation works with the latest commit of BOUT++ version 3

The following folders suffers from the `0.51 bug` (the perpendicular
resistivity multiplied with `0.51` rather than `1.0`).

* additionalPlots: Additional plots related to the CELMA code
* celma: Folder with the codes - Suffers from `0.51 bug`
* celmaCurMom: CELMA codes where the current and parallel momentum is solved -
  Suffers from `0.51 bug`
* celmaApar: CELMA codes where the $A_\|$ is accounted for -
  Suffers from `0.51 bug`
* celmaRC: CELMA release candidate
* derivation: Derivation of operators and boundaries
* MES: Checking if implementation is convergent
* parameters: Calculation of typical parameters
* toolbox: Scripts and functions which are not essential to the runs

Branch: master - Running and stable with latest BOUT-dev commit (BOUT++ 3)
