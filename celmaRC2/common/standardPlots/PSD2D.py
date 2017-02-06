#!/usr/bin/env python

"""Post-processor for the PSD2D"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.PSD import DriverPSD

#{{{PSD2DPlot
def PSD2DPlot(dmp_folders, collectPaths, plotSuperKwargs, tSlice = None):
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

    varName = "n"
    mode    = "fluct"

    yInd    = 16
    zInd    = 0
    tSlice  = tSlice

    indicesArgs   = (None, yInd, zInd)
    indicesKwargs = {"tSlice" : tSlice}

    plotLimits = {"xlim":None      ,\
                  "ylim":(100, 3e4),\
                  "zlim":(-7,0)}

    dPSD = DriverPSD(
                     # DriverPSD
                     dmp_folders                   ,\
                     indicesArgs                   ,\
                     indicesKwargs                 ,\
                     plotSuperKwargs               ,\
                     varName           = varName   ,\
                     mode              = mode      ,\
                     plotLimits        = plotLimits,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dPSD.driverPSD2D()
#}}}
