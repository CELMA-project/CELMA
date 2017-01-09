#!/usr/bin/env python

"""Post-processor for performance"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.performance import DriverPerformance

#{{{performancePlot
def performancePlot(dmp_folders, collectPaths, plotSuperKwargs,\
                    allFolders = False):
    #{{{docstring
    """
    Runs the standard performance plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    allFolders : bool
        If "init" and "expand" has been included in the plot.
    """
    #}}}

    useSubProcess = False
    convertToPhysical = True

    dP = DriverPerformance(
                     # DriverPerformance
                     dmp_folders            ,\
                     convertToPhysical      ,\
                     plotSuperKwargs        ,\
                     allFolders = allFolders,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dP.driverPerformance()
#}}}
