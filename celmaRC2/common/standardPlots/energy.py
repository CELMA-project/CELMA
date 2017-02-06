#!/usr/bin/env python

"""Post-processor for energy"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.energy import DriverEnergy

#{{{energyPlot
def energyPlot(dmp_folders, collectPaths, plotSuperKwargs, tSlice = None):
    #{{{docstring
    """
    Runs the standard energy plot

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

    if tSlice is not None:
        sliced = True
    else:
        sliced = False

    plotSuperKwargs["sliced"] = sliced

    dE = DriverEnergy(
                     # DriverEnergy
                     dmp_folders    ,\
                     tSlice         ,\
                     plotSuperKwargs,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dE.driverEnergy()
#}}}
