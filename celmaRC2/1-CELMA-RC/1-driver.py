#!/usr/bin/env python

"""Driver which checks the plots."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.scanDriver import ScanDriver
from CELMAPy.driverHelpers import pathMerger
from bout_runners import basic_runner
from standardPlots import fields1DAnimation
import pickle

# DRIVE
# =============================================================================
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

scanB0.setInitOptions(timestep = 1e-10, nout = 2)

scanB0.setRunTypeOptions(runInit   = True ,\
                         runExpand = False,\
                         runLin    = False,\
                         runTurb   = False)



# Set common runner options
scanB0.setCommonRunnerOptions(nproc      = 4   ,\
                              cpy_source = True,\
                             )

# Run
scanB0.runScan()
# =============================================================================


# POST-PROCESS
# =============================================================================
# Have support $\parallel$
import matplotlib.pylab as plt
plt.rc("text", usetex=False)
with open(os.path.join(directory, "dmpFoldersDict.pickle"), "rb") as f:
        dmpFolders = pickle.load(f)

collectPaths = pathMerger(dmpFolders, ("init",))["param0"]
plotSuperKwargs = {
                    "showPlot"        : False        ,\
                    "savePlot"        : True         ,\
                    "savePath"        : None         ,\
                    "savePathFunc"    : "onlyScan"   ,\
                    "extension"       : None         ,\
                    "timeStampFolder" : False        ,\
                    # scanParameter needed in onlyScans
                    "scanParameter"   : "B0"         ,\
                   }

fields1DAnimation((dmpFolders["init"][0],), collectPaths, plotSuperKwargs)
# =============================================================================
