#!/usr/bin/env python

"""Post-processor for energy"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.energy import DriverEnergy

#{{{energyPlot
def energyPlot(dmp_folders, collectPaths, tSlice = None):
    #{{{docstring
    """
    Runs the standard energy plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useSubProcess     = False
    convertToPhysical = True

    if tSlice is not None:
        sliced = True
    else:
        sliced = False

    plotSuperKwargs = {\
                        "showPlot"        : False ,\
                        "savePlot"        : True  ,\
                        "savePath"        : None  ,\
                        "savePathFunc"    : None  ,\
                        "extension"       : None  ,\
                        "dmp_folders"     : None  ,\
                        "timeStampFolder" : False ,\
                        "sliced"          : sliced,\
                       }

    dE = DriverEnergy(
                     # DriverEnergy
                     dmp_folders    ,\
                     tSlice         ,\
                     plotSuperKwargs,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dE.driverEnergy()
#}}}
