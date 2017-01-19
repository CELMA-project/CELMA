#!/usr/bin/env python

"""Post-processor for the total flux"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.totalFlux import DriverTotalFlux

#{{{totalFluxPlot
def totalFluxPlot(dmp_folders, collectPaths, plotSuperKwargs, tSlice = None):
    #{{{docstring
    """
    Runs the standard total flux plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders.
    collectPaths : tuple
        Tuple of the paths to collect from.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useSubProcess     = False
    xInd              = None
    yInd              = None
    mode              = "normal"
    convertToPhysical = True

    dTF = DriverTotalFlux(\
                     # DriverTotalFlux
                     dmp_folders                          ,\
                     plotSuperKwargs                      ,\
                     xInd              = xInd             ,\
                     yInd              = yInd             ,\
                     tSlice            = tSlice           ,\
                     mode              = mode             ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dTF.driverTotalFlux()
#}}}
