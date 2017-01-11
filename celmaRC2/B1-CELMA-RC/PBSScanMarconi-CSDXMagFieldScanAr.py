#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.scanDriver import ScanDriver

directory = "CSDXMagFieldScanAr"

# Create object
scanB0 = ScanDriver(directory)

# Set the scan
# NOTE: The scan must be in descending order in order for the growth
#       rates post-processing to work
B0 = (  1.0e-1,   8.0e-2,   6.0e-2,   4.0e-2,  2.0e-2)
Lx = (  7.8633,   6.2906,   4.7180,   3.1453,  1.5727)
Ly = (275.2144, 220.1715, 165.1286, 110.0858, 55.0429)

scanParameters  = ("B0", "Lx", "Ly")
series_add = (\
              ("input", "B0", B0),\
              ("geom" , "Lx", Lx),\
              ("geom" , "Ly", Ly),\
             )

# Set the main options
scanB0.setMainOptions(\
                       scanParameters   = scanParameters,\
                       series_add       = series_add    ,\
                       theRunName       = directory     ,\
                       make             = False         ,\
                       boutRunnersNoise = {"vort":1e-6} ,\
                     )

scanB0.setInitOptions(BOUT_walltime = "72:00:00")

scanB0.setExpandOptions(timestep      = 25       ,\
                        nout          = 2        ,\
                        BOUT_walltime = "72:00:00")

# Set common runner options
scanB0.setCommonRunnerOptions(\
                              nproc        = 48            ,\
                              cpy_source   = True          ,\
                              BOUT_nodes   = 2             ,\
                              BOUT_ppn     = 36            ,\
                              BOUT_queue   = "xfualongprod",\
                              BOUT_account = "FUA11_SOLF"  ,\
                             )

# Run
scanB0.runScan(restartTurb=3, boussinesq = True)
