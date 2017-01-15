#!/usr/bin/env python

"""Driver which checks the plots."""

from boututils.options import BOUTOptions
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
directory = "BousCSDXNeutralScanAr"

# Create object
scanNn = ScanDriver(directory, runner = basic_runner)

# Set the scan
# NOTE: The scan must be in descending order in order for the growth
#       rates post-processing to work
ionizationDegrees = (1,)
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

# Set the options
scanNn.setMainOptions(\
                       scanParameters   = scanParameters,\
                       series_add       = series_add    ,\
                       theRunName       = directory     ,\
                       make             = True          ,\
                       boutRunnersNoise = {"vort":1e-6} ,\
                     )

scanNn.setInitOptions(timestep = 1e-10, nout = 2)

scanNn.setRunTypeOptions(runInit   = True ,\
                         runExpand = False,\
                         runLin    = False,\
                         runTurb   = False)



# Set common runner options
scanNn.setCommonRunnerOptions(nproc      = 4   ,\
                              cpy_source = True,\
                             )

# Run
scanNn.runScan(boussinesq = True)
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
                    "scanParameter"   : "nn"         ,\
                   }

fields1DAnimation((dmpFolders["init"][0],), collectPaths,\
        plotSuperKwargs, boussinesq=True)
# =============================================================================
