#!/usr/bin/env python

"""
Example which fixes the restart files after premature exit.
To be run after restartRunsExample.py

**NOTE**: DO NOT RUN THIS unless tere is a problem.
"""

import glob
from shutil import copy2

# INPUT
# =============================================================================
# Set the scan
B0s = (  1.0e-1,   8.0e-2,   6.0e-2)
Lxs = (  7.8633,   6.2906,   4.7180)
Lys = (275.2144, 220.1715, 165.1286)
# New is here referring to the restarted runs
newNout     = 1000
newTimestep = 2.0
# Old is here referring to the runs where the premature exit occured
oldNout     = 2
oldTimestep = 2000.0
# =============================================================================

for B0, Lx, Ly in zip(B0s, Lxs, Lys):
    pathFrom = "CSDXMagFieldScanAr/nout_{}_timestep_{}/nz_1/geom_Lx_{}_geom_Ly_{}_input_B0_{}_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/".\
            format(newNout, newTimestep, Lx, Ly, B0)
    pathTo = "CSDXMagFieldScanAr/nout_{}_timestep_{}/nz_1/geom_Lx_{}_geom_Ly_{}_input_B0_{}_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/".\
            format(oldNout, oldTimestep, Lx, Ly, B0)

    restartFiles = glob.glob(pathFrom + "*.restart.*")

    for f in restartFiles:
        copy2(f, pathTo)
