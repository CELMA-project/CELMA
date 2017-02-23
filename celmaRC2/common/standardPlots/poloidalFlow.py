#!/usr/bin/env python

"""Post-processor for the poloidal flow"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.poloidalFlow import DriverPoloidalFlow

#{{{poloidalFlowPlot
def poloidalFlowPlot(dmp_folders, collectPaths, steadyStatePath, plotSuperKwargs,\
                  tSlice = None):
    #{{{docstring
    """
    Runs the standard 2d spectral density plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useMultiProcess     = False
    convertToPhysical = True
    yInd              = 16
    tSlice            = tSlice

    dRP = DriverPoloidalFlow(
                     # DriverPoloidalFlow
                     dmp_folders                          ,\
                     steadyStatePath                      ,\
                     yInd                                 ,\
                     tSlice                               ,\
                     plotSuperKwargs                      ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dRP.driverPoloidalFlow()
#}}}
