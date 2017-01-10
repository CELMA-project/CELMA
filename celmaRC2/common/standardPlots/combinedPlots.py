#!/usr/bin/env python

"""Post-processor for combinedPlots"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.combinedPlots import DriverCombinedPlots

#{{{combinedPlotsPlot
def combinedPlotsPlot(dmp_folders, collectPaths, steadyStatePath,\
                      plotSuperKwargs, tSlice = None):
    #{{{docstring
    """
    Runs the standard combined plots plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    steadyStatePath : str
        The corresponding steady state path
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useSubProcess = False

    varName           = "n"
    convertToPhysical = True
    mode              = "fluct"

    yInd              = 16
    zInd              = 0

    if tSlice is not None:
        sliced = True
    else:
        sliced = False

    useSubProcess = False

    dTT = DriverCombinedPlots(
                     # DriverCombinedPlots
                     dmp_folders                ,\
                     steadyStatePath            ,\
                     yInd                       ,\
                     zInd                       ,\
                     tSlice                     ,\
                     plotSuperKwargs            ,\
                     varName           = varName,\
                     mode              = mode   ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dTT.driverCombinedPlots()
#}}}
