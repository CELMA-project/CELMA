#!/usr/bin/env python

"""Post-processor for blobs"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.blobs import DriverBlobs

UNDER CONSTRUCTION
#{{{blobsPlot
def blobsPlot(dmp_folders,\
              collectPaths,\
              steadyStatePath,\
              plotSuperKwargs,\
              tSlice = None):
    #{{{docstring
    """
    Runs the standard blobs plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders.
    collectPaths : tuple
        Tuple of the paths to collect from.
    steadyStatePath : str
        The corresponding steady state path.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : [None|Slice]
        How to slice the time.
    """
    #}}}

    useSubProcess = False

    varName           = "n"
    convertToPhysical = True

    xInd              = None
    yInd              = 16
    zInd              = None
    nPoints           = 3
    equallySpace      = "x"

    indicesArgs   = (xInd, yInd, zInd)
    indicesKwargs = {"tSlice"          : tSlice         ,\
                     "nPoints"         : nPoints        ,\
                     "equallySpace"    : equallySpace   ,\
                     "steadyStatePath" : steadyStatePath,\
                     }

    dMS = DriverBlobs(
                     # DriverBlobs
                     dmp_folders              ,\
                     indicesArgs              ,\
                     indicesKwargs            ,\
                     plotSuperKwargs          ,\
                     varName         = varName,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dMS.driverBlobs()
#}}}
