#!/usr/bin/env python

"""Driver which runs using PBS."""

from boututils.options import BOUTOptions
import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.scanDriver import ScanDriver

directory = "CSDXNeutralScanAr"

# Create object
scanNn = ScanDriver(directory)

# Set the scan
# NOTE: The scan must be in descending order in order for the growth
#       rates post-processing to work
ionizationPercents = (80, 60, 40, 20)
option = BOUTOptions(directory)
n0 = float(option.input["n0"])
nn = tuple(n0/((pct)/100) for pct in ionizationPercents)

scanParameters  = ("nn",)
series_add = (\
              ("input", "nn", nn),\
             )

# Set the main options
scanNn.setMainOptions(\
                       scanParameters   = scanParameters,\
                       series_add       = series_add    ,\
                       theRunName       = directory     ,\
                       make             = False         ,\
                       boutRunnersNoise = {"vortD":1e-6},\
                     )

scanNn.setInitOptions(BOUT_walltime = "72:00:00")

scanNn.setExpandOptions(timestep      = 25       ,\
                        nout          = 2        ,\
                        BOUT_walltime = "72:00:00")

# Set common runner options
scanNn.setCommonRunnerOptions(\
                              nproc        = 48            ,\
                              cpy_source   = True          ,\
                              BOUT_nodes   = 2             ,\
                              BOUT_ppn     = 36            ,\
                              BOUT_queue   = "xfualongprod",\
                              BOUT_account = "FUA11_SOLF"  ,\
                             )

# Run
scanNn.runScan(restartTurb=3)
