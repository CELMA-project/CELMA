#!/usr/bin/env python

"""Post-processor for fourierModes"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fourierModes import DriverFourierModes

#{{{standardFourierModesPlot
def standardFourierModesPlot(dmp_folders, collectPaths, tSlice = None):
    """
    Runs the standard fourier modes plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    tSlice : [None|Slice]
        How to slice the time.
    """

    useSubProcess = True

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd              = None
    yInd              = 16
    zInd              = None
    nPoints           = 1
    equallySpace      = "x"

    steadyStatePath   = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    indicesArgs   = (xInd, yInd, zInd)
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

    print("\n\nTesting fourier mode driver")
    dFM = DriverFourierModes(
                     # DriverFourierModes
                     dmp_folders                ,\
                     indicesArgs                ,\
                     indicesKwargs              ,\
                     plotSuperKwargs            ,\
                     varName           = varName,\
                     nModes            = nModes ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dFM.driverFourierMode()
    print("Success!\n\n")
#}}}
