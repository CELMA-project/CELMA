# CELMA-dev

Developement branch of the CELMA code, copied from revision 1130 of
https://svn.fysik.dtu.dk/mmag/Code-developement/BOUT/

Current implementation works with the latest commit of BOUT++ version 3

**NOTE**:
1. The `0.51 bug` - the perpendicular resistivity multiplied with `0.51` rather
   than `1.0`.
2. The `parallel bug` - Mistake in momentum density derivation spreads to
   current and sum of momentum density equation in the code.

* additionalPlots: Additional plots related to the CELMA code
* celma: Folder with the codes - Suffers from `0.51 bug` and `parallel bug`
* celmaCurMom: CELMA codes where the current and parallel momentum is solved -
  Suffers from `0.51 bug` and `parallel bug`
* celmaApar: CELMA codes where the $A_\|$ is accounted for -
  Suffers from `0.51 bug` and `parallel bug`
* celmaRC: CELMA release candidate - Suffers from `parallel bug`
* celmaRC2: CELMA release candidate 2
* derivation: Derivation of operators and boundaries
* MES: Checking if implementation is convergent
* parameters: Calculation of typical parameters
* toolbox: Scripts and functions which are not essential to the runs

Branch: master - Running and stable with latest BOUT-dev commit (BOUT++ 3)
