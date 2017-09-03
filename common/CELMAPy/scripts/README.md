# scripts

Contains non-callable scripts.
To run the scripts, copy the `.py` files to the root of the project and
run the script.

Example:

```bash
~/CELMA/common/CELMAPy/scripts $ scp prematureExitFixes/restartRunsExample.py ../../../celma/
~/CELMA/common/CELMAPy/scripts $ cd ../../../celma/
~/CELMA/celma $ python restartRunsExample.py
```

[prematureExitFixes](prematureExitFixes) - Contains examples of how to fix premature exits
[captureAllRestartAndLogFiles.py](captureAllRestartAndLogFiles.py) - Saves all *.log.* and *.restart.* files of a directory to a zip.
[refreshDates.py](refreshDates.py) - Refresh dates of files to prevent automatic deletion  by cluster.
[refreshDatesPBSDriver.py](refreshDatesPBSDriver.py) - Submits refreshDates() to the PBS queue.
