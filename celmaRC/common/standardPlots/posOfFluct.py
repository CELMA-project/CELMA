#!/usr/bin/env python

"""Post-processor for the position of the radial profile"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.radialProfile import DriverRadialProfile

#{{{posOfFluctPlot
def posOfFluctPlot(dmp_folders, collectPaths, steadyStatePath, tSlice = None):
    #{{{docstring
    """
    Runs the standard position of the radial profile plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    steadyStatePath : str
        Path to the steady state
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useSubProcess     = False
    convertToPhysical = True

    plotSuperKwargs = {\
                        "showPlot"        : False ,\
                        "savePlot"        : True  ,\
                        "savePath"        : None  ,\
                        "savePathFunc"    : None  ,\
                        "extension"       : None  ,\
                        "dmp_folders"     : None  ,\
                        "timeStampFolder" : False ,\
                       }

    varName  = "n"
    var2Name = "phi"
    yInd     = 16
    tSlice   = tSlice

    dRP = DriverRadialProfile(
                     # DriverRadialProfile
                     dmp_folders                          ,\
                     steadyStatePath                      ,\
                     yInd                                 ,\
                     tSlice                               ,\
                     plotSuperKwargs                      ,\
                     varName           = varName          ,\
                     var2Name          = var2Name         ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dRP.driverPosOfFluct()
#}}}
