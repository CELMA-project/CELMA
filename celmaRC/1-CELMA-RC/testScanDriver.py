#!/usr/bin/env python

"""Test of the ScanDriver class"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.scanDriver import ScanDriver

directory = "CSDXMagFieldScan"

# Create object
# NOTE: This is a test using basic_runner
from bout_runners import basic_runner
scanB0 = ScanDriver(directory, runner = basic_runner)

# Set the scan
B0 = (   1.0e-1,    8.0e-2,    6.0e-2,   4.0e-2,   2.0e-2)
Lx = (  49.5195,   39.6156,   29.7117,  19.8078,   9.9039)
Ly = (1733.1827, 1386.5461, 1039.9096, 693.2731, 346.6365)

scanParameters  = ("B0", "Lx", "Ly")
series_add = (\
              ("input", "B0", B0),\
              ("geom" , "Lx", Lx),\
              ("geom" , "Ly", Ly),\
             )


# Set the main options
scanB0.setMainOptions(\
                       scanParameters = scanParameters  ,\
                       series_add     = series_add      ,\
                       theRunName     = directory       ,\
                       make           = False           ,\
                       boutRunnersNoise = {"vortD":1e-6},\
                     )

# Set common runner options
scanB0.setCommonRunnerOptions(\
                              nproc      = 48  ,\
                              cpy_source = True,\
                              BOUT_nodes = 3   ,\
                              BOUT_ppn   = 16  ,\
                             )

# Run
scanB0.runScan(restartTurb=2)
