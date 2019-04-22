#!/usr/bin/env python

"""
Example script which fixes the log files after premature exit.
To be run after restartRunsExample.py

**NOTE**: DO NOT RUN THIS unless tere is a problem.
"""

import glob, datetime
import numpy as np
import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.logReader import getLogNumbers

# INPUT
# =============================================================================
# How many more simtime entries should there have been in the log files
missingEntries = 2
# New is here referring to the restarted runs
newNout     = 1000
newTimestep = 2.0
# Old is here referring to the runs where the premature exit occured
oldNout     = 2
oldTimestep = 2000.0
# Set the scan
B0s = (  1.0e-1,   8.0e-2,   6.0e-2)
Lxs = (  7.8633,   6.2906,   4.7180)
Lys = (275.2144, 220.1715, 165.1286)
# =============================================================================

for B0, Lx, Ly in zip(B0s, Lxs, Lys):
    pathFrom = "CSDXMagFieldScanAr/nout_{}_timestep_{}/nz_1/geom_Lx_{}_geom_Ly_{}_input_B0_{}_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/".\
            format(newNout, newTimestep, Lx, Ly, B0)
    pathTo = "CSDXMagFieldScanAr/nout_{}_timestep_{}/nz_1/geom_Lx_{}_geom_Ly_{}_input_B0_{}_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/".\
            format(oldNout, oldTimestep, Lx, Ly, B0)

    log = getLogNumbers(pathFrom)

    tmps = [None]*missingEntries
    for key in log.keys():
        # Create the desired number of entries
        tmps = np.array_split(log[key], missingEntries)
        # Reset the log
        log[key] = []
        if ("Time" not in key) and ("RHSevals" not in key):
            for t in tmps:
                log[key].append(t.mean())
        elif key == "SimTime":
            for t in tmps:
                # Only the last simulation time counts
                log[key].append(t[-1])
        elif (key == "WallTime") or (key == "RHSevals"):
            # Walltime and RHS evals are cumulative
            for t in tmps:
                log[key].append(t.sum())

    logFiles = glob.glob(pathTo + "*.log.*")

    for lf in logFiles:
        # Print the new entries
        with open(lf, "a") as f:
            for e in range(missingEntries):
                string =\
                    ("{{SimTime[{0}]:.3e}}   {{RHSevals[{0}]:.0f}}   "
                     "{{WallTime[{0}]:.2e}}   {{Calc[{0}]:.1f}}   "
                     "{{Inv[{0}]:.1f}}   {{Comm[{0}]:.1f}}   {{I/O[{0}]:.1f}}   "
                     "{{SOLVER[{0}]:.1f}}\n").format(e)
                f.write(string.format(**log))

        # Read the old Walltime
        with open(lf, "r") as f:
            oldLog = getLogNumbers(pathTo)
            oldWallTime = sum(oldLog["WallTime"])

        # Not entirely the wall time, as the oldWallTime represents a mean
        # of all the files
        approxWallTime = oldWallTime + sum(log["WallTime"])

        with open(lf, "a") as f:
            time = datetime.timedelta(seconds = approxWallTime)
            f.write("\n!!!This log file has been patched!!!\n")
            f.write("\nRun time approximately {}\n".format(str(time)))
