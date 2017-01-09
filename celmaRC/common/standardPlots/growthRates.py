#!/usr/bin/env python

"""Post-processor for growthRates"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.growthRates import DriverGrowthRates

#{{{growthRatesPlot
def growthRatesPlot(dmp_folders     ,\
                    scanCollectPaths,\
                    steadyStatePaths,\
                    scanParameter   ,\
                    tSlices         ,\
                    plotSuperKwargs ,\
                   ):
    #{{{docstring
    """
    Runs the standard growth rates plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple of tuple
        Tuple of the paths to collect from (one for each scan)
    steadyStatePaths : tuple
        The corresponding steady state paths
    scanParameter : str
        Name of the scan parameter (as in the scan driver)
    tSlices : tuple of slices
        The time slices to use for each folder in the scan order.
        This must manually be selected to the linear phase found in
        from the fourier moedes.
        The tuple must be ordered according to the scan order in
        scanCollectPaths.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd            = None
    yInd            = 16
    nPoints         = 3
    equallySpace    = "x"
    steadyStatePath = None

    useSubProcess = True

    indicesArgs   = (xInd, yInd)
    indicesKwargs = {"tSlice"          : None           ,\
                     "nPoints"         : nPoints        ,\
                     "equallySpace"    : equallySpace   ,\
                     "steadyStatePath" : steadyStatePath,\
                     }

    dGR = DriverGrowthRates(
                     # DriverGrowthRates
                     dmp_folders                          ,\
                     scanCollectPaths                     ,\
                     steadyStatePaths                     ,\
                     tSlices                              ,\
                     scanParameter                        ,\
                     indicesArgs                          ,\
                     indicesKwargs                        ,\
                     plotSuperKwargs                      ,\
                     varName           = varName          ,\
                     nModes            = nModes           ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     useSubProcess = useSubProcess,\
                          )
    dGR.driverGrowthRates()
#}}}
