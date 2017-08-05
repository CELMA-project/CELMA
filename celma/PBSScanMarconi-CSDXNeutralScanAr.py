#!/usr/bin/env python

"""Driver which runs using PBS."""

from boututils.options import BOUTOptions
import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.scanDriver import ScanDriver

# INPUT
# =============================================================================
# If the queuing system uses accounts, set the account here
# Example: "FUA11_SOLF"
account = None
# Usually, the queueing system has its own default queue, if not,
# specify here
# Example: "xfualongprod"
queue = None
# If you would like a mail on finished job enter your mail here
# Example: "john@doe.com"
mail = None
# =============================================================================

directory = "CSDXNeutralScanAr"

# Create object
scanNn = ScanDriver(directory)

# Set the scan
# NOTE: The scan must be in descending order in order for the growth
#       rates post-processing to work
ionizationDegrees = (80, 60, 40, 20, 1)
option = BOUTOptions(directory)
n0 = float(option.input["n0"])
# Ionization degree
# d = ni/(ni+nn) => nn = (ni/d) - ni
# NOTE: Conversion from percent to fraction
nn = tuple(n0/(d/100) - n0 for d in ionizationDegrees)

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
                              nproc        = 48     ,\
                              cpy_source   = True   ,\
                              BOUT_nodes   = 2      ,\
                              BOUT_ppn     = 36     ,\
                              BOUT_queue   = queue  ,\
                              BOUT_account = account,\
                              BOUT_mail    = mail   ,\
                             )

# Run
scanNn.runScan(restartTurb=2)
