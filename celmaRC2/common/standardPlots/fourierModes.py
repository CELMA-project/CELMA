#!/usr/bin/env python

"""Post-processor for fourierModes"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fourierModes import DriverFourierModes

#{{{fourierModesPlot
def fourierModesPlot(dmp_folders,\
                     collectPaths,\
                     steadyStatePath,\
                     plotSuperKwargs,\
                     tSlice = None):
    #{{{docstring
    """
    Runs the standard fourier modes plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders.
    collectPaths : tuple
        Tuple of the paths to collect from.
    steadyStatePath : str
        The corresponding steady state path.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useSubProcess = False

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd              = None
    yInd              = 16
    zInd              = None
    nPoints           = 1
    equallySpace      = "x"

    if tSlice is not None:
        sliced = True
    else:
        sliced = False

    plotSuperKwargs["sliced"] = sliced

    indicesArgs   = (xInd, yInd, zInd)
    indicesKwargs = {"tSlice"          : tSlice         ,\
                     "nPoints"         : nPoints        ,\
                     "equallySpace"    : equallySpace   ,\
                     "steadyStatePath" : steadyStatePath,\
                     }

    dFM = DriverFourierModes(
                     # DriverFourierModes
                     dmp_folders                ,\
                     indicesArgs                ,\
                     indicesKwargs              ,\
                     plotSuperKwargs            ,\
                     varName           = varName,\
                     nModes            = nModes ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dFM.driverFourierMode()
#}}}
