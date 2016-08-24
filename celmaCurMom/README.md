# celmaCurMom

Folder with the different versions of the code when evolving the parallel current and
total parallel density momentum.

Each run of each group has its own driver. `bout_runners` is used for running
and post-processing.

* 1-CELMACurMom - Original code, copied from `celma/8.3-CELMASplitCleanUp`
* 2-CELMACurMomWMonitors - As 1, but added monitors
* common - python post processing and own implementations to BOUT++
