#!/usr/bin/env python

"""Post-processor for performance"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.performance import DriverPerformance

#{{{performancePlot
def performancePlot(dmp_folders, collectPaths, allFolders = False):
    #{{{docstring
    """
    Runs the standard performance plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    allFolders : bool
        If "init" and "expand" has been included in the plot.
    """
    #}}}

    useSubProcess = False
    convertToPhysical = True

    plotSuperKwargs = {\
                        "showPlot"        : False     ,\
                        "savePlot"        : True      ,\
                        "savePath"        : None      ,\
                        "extension"       : None      ,\
                        # NOTE: No implemented func which doesn't
                        #       require theRunName yet
                        "savePathFunc"    : None      ,\
                        "dmp_folders"     : None      ,\
                        "timeStampFolder" : False     ,\
                        "allFolders"      : allFolders,\
                       }

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
