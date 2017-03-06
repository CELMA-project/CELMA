#!/usr/bin/env python

"""Post-processor for performance"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.performance import DriverPerformance

#{{{performancePlot
def performancePlot(dmp_folders, collectPaths, mode, plotSuperKwargs, tSlice=None):
    #{{{docstring
    """
    Runs the standard performance plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    mode : ["init"|"expand"|"linear"|"turbulence"|"all"]
        What part of the simulation is being plotted for.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : slice
        Use if the data should be sliced.
    """
    #}}}

    useMultiProcess = False
    convertToPhysical = True

    dP = DriverPerformance(
                     # DriverPerformance
                     dmp_folders      ,\
                     convertToPhysical,\
                     mode             ,\
                     plotSuperKwargs  ,\
                     tSlice = tSlice  ,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dP.driverPerformance()
#}}}
