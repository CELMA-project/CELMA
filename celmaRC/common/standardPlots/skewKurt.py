#!/usr/bin/env python

"""Post-processor for the skewness and kurtosis"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.skewnessKurtosis import DriverSkewnessKurtosis

#{{{skewKurtPlot
def skewKurtPlot(dmp_folders, collectPaths, tSlice = None):
    #{{{docstring
    """
    Runs the standard 2d spectral density plot

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

    varName           = "n"
    mode              = "fluct"

    yInd              = 16
    zInd              = 0
    tSlice            = None

    indicesArgs   = (None, yInd, zInd)
    indicesKwargs = {"tSlice" : tSlice}

    plotSuperKwargs = {\
                        "showPlot"        : False,\
                        "savePlot"        : True ,\
                        "savePath"        : None ,\
                        "savePathFunc"    : None ,\
                        "extension"       : None ,\
                        "dmp_folders"     : None ,\
                        "timeStampFolder" : False,\
                       }

    dSK = DriverSkewnessKurtosis(
                     # DriverSkewnessKurtosiss
                     dmp_folders                   ,\
                     indicesArgs                   ,\
                     indicesKwargs                 ,\
                     plotSuperKwargs               ,\
                     varName           = varName   ,\
                     mode              = mode      ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dSK.driverSkewnessKurtosis()
#}}}
