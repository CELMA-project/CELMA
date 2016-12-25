#!/usr/bin/env python

"""Post-processor test for timeTraces"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.timeTrace import driverTimeTrace

#{{{timeTraceTest
def timeTraceTest():
    """
    Runs the single perp test for the time traces
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    varName           = "lnN"
    convertToPhysical = True
    mode              = "fluct"

    xInd              = 16
    yInd              = 16
    zInd              = 128
    tSlice            = None
    nPoints           = 3
    equallySpace      = None
    steadyStatePath   = None

    indicesArgs   = (xInd, yInd, zInd)
    indicesKwargs = {"tSlice"          : tSlice         ,\
                     "nPoints"         : nPoints        ,\
                     "equallySpace"    : equallySpace   ,\
                     "steadyStatePath" : steadyStatePath,\
                     }

    savePath          = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting single time trace")
    driverTimeTrace(collectPaths     ,\
                    varName          ,\
                    convertToPhysical,\
                    mode             ,\
                    indicesArgs      ,\
                    indicesKwargs    ,\
                    plotSuperKwargs  ,\
                   )

    print("Success!\n\n")
#}}}



#{{{driverTest
def driverTest():
    """
    Runs the driver test for the time traces
    """
    dmp_folders  = (
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0",
            )
    collectPaths =\
        (\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0",
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",
        )

# FIXME: Try n
    varName           = "lnN"
    convertToPhysical = True
    savePathFunc    = "scanWTagSaveFunc"
    useSubProcess     = False
    theRunName        = "test"
    mode              = "fluct"
    nPoints           = 5
    xInd              = None
    yInd              = 16
    zInd              = 0
    tSlice            = None

    equallySpace = "x"
    steadyStatePath   = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    print("\n\nTesting time trace driver")
    dtt = DriverTimeTrace(
                # DriverPostProcessingSuperClass
                dmp_folders                          ,\
                collectPaths      = collectPaths     ,\
                convertToPhysical = convertToPhysical,\
                showPlot          = False            ,\
                savePlot          = True             ,\
                saveFolder        = None             ,\
                savePathFunc      = savePathFunc     ,\
                useSubProcess     = useSubProcess    ,\
                extension         = "png"            ,\
                # scanWTagSaveFunc
                theRunName = theRunName,\
                # CollectAndCalcPointsSuperClass
                mode            = mode           ,\
                nPoints         = nPoints        ,\
            )
    dtt.setVarName(varName)
    dtt.setIndices(xInd, yInd, zInd,\
                   tSlice = tSlice, equallySpace = equallySpace,\
                   steadyStatePath = steadyStatePath)
    import pdb; pdb.set_trace()
    dtt.plotTimeTrace()
    print("\n\nSuccess!")
#}}}

if __name__ == "__main__":
    driverTest()
