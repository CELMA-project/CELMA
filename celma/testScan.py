#!/usr/bin/env python

"""Driver which checks that celma is properly running."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.scanDriver import ScanDriver
from bout_runners import basic_runner

directory = "CSDXMagFieldScanAr"

# Create object
scanB0 = ScanDriver(directory, runner = basic_runner)

# Set the scan
# NOTE: The scan must be in descending order in order for the growth
#       rates post-processing to work
B0 = (  1.0e-1)
Lx = (  7.8633)
Ly = (275.2144)

scanParameters  = ("B0", "Lx", "Ly")
series_add = (\
              ("input", "B0", B0),\
              ("geom" , "Lx", Lx),\
              ("geom" , "Ly", Ly),\
             )

# Set the options
scanB0.setMainOptions(\
                       scanParameters   = scanParameters,\
                       series_add       = series_add    ,\
                       theRunName       = directory     ,\
                       make             = False         ,\
                       boutRunnersNoise = {"vortD":1e-6},\
                     )

scanB0.setInitOptions      (timestep = 1e-10, nout = 2)
scanB0.setExpandOptions    (timestep = 1e-10, nout = 2)
scanB0.setLinearOptions    (timestep = 1e-10, nout = 2)
scanB0.setTurbulenceOptions(timestep = 1e-10, nout = 2)

# Set common runner options
scanB0.setCommonRunnerOptions(nproc = 4, cpy_source = True)

# Run
scanB0.runScan()
