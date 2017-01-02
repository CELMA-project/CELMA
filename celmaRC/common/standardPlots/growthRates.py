#!/usr/bin/env python

"""Post-processor for growthRates"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.growthRates import DriverGrowthRates

#{{{growthRatesPlot
def growthRatesPlot(dmp_folders, scanCollectPaths, steadyStatePaths, scanParameter, startInds, endInds):
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
    startInds : tuple
        Tuple of the start indices
    startInds : tuple
        Tuple of the end indices
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    scanParameter = "B0"

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd            = None
    yInd            = 16
    tSlice          = None
    nPoints         = 3
    equallySpace    = "x"
    steadyStatePath = None

    useSubProcess = True

    indicesArgs   = (xInd, yInd)
    indicesKwargs = {"tSlice"          : tSlice         ,\
                     "nPoints"         : nPoints        ,\
                     "equallySpace"    : equallySpace   ,\
                     "steadyStatePath" : steadyStatePath,\
                     }

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    dGR = DriverGrowthRates(
                     # DriverGrowthRates
                     dmp_folders                          ,\
                     scanCollectPaths                     ,\
                     steadyStatePaths                     ,\
                     startInds                            ,\
                     endInds                              ,\
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
