#!/usr/bin/env python

"""Post-processor for the blob density PDF"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.PDF import DriverPDF
from .blobs import xInd, yInd, zInd

#{{{blobDensPDF
def blobDensPDF(dmp_folders, collectPaths, plotSuperKwargs, tSlice = None):
    #{{{docstring
    """
    Runs the standard blob density PDF.

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

    useMultiProcess   = False
    convertToPhysical = True
    varName           = "n"
    mode              = "fluct"

    indicesArgs   = (xInd, yInd, zInd)
    nPoints       = 1
    indicesKwargs = {"tSlice" : tSlice, "nPoints" : nPoints}

    dPDF =  DriverPDF(
                     # DriverPDFs
                     dmp_folders                ,\
                     indicesArgs                ,\
                     indicesKwargs              ,\
                     plotSuperKwargs            ,\
                     varName           = varName,\
                     mode              = mode   ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dPDF.driverPDF()
#}}}
