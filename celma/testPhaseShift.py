#!/usr/bin/env python

"""Post-processor test for phaseShift"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.growthRates import DriverGrowthRates, DriverPhaseShift, driverPhaseShift

#{{{phaseShiftTest
def phaseShiftTest():
    """
    Runs the test for the growth rates
    """

    scanCollectPaths =\
       (\
        (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/"\
        ),\
        (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
        ),\
       )

    steadyStatePaths =\
        (\
         "CSDXMagFieldScanAr/nout_2_timestep_25/nz_256/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
         "CSDXMagFieldScanAr/nout_2_timestep_25/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
        )

    tSlices = (slice(80,240), slice(80,210))

    scanParameter = "B0"

    collectArgs = DriverGrowthRates.makeCollectArgs(scanCollectPaths,\
                                                    steadyStatePaths,\
                                                    tSlices         ,\
                                                    scanParameter)

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd            = 16
    yInd            = 16
    zInd            = 0
    tSlice          = None
    nPoints         = 1
    equallySpace    = "x"
    steadyStatePath = None

    indicesArgs   = (xInd, yInd, zInd)
    indicesKwargs = {"tSlice"          : tSlice         ,\
                     "nPoints"         : nPoints        ,\
                     "equallySpace"    : equallySpace   ,\
                     "steadyStatePath" : steadyStatePath,\
                     }

    getDataArgs = DriverGrowthRates.makeGetDataArgs(varName          ,\
                                                    convertToPhysical,\
                                                    indicesArgs      ,\
                                                    indicesKwargs    ,\
                                                    nModes)

    savePath          = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting phase shifts")
    driverPhaseShift(collectArgs    ,\
                     getDataArgs    ,\
                     plotSuperKwargs,\
                    )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the growth rates
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)

    scanCollectPaths =\
       (\
        (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/"\
        ),\
        (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
        ),\
       )

    steadyStatePaths =\
        (\
         "CSDXMagFieldScanAr/nout_2_timestep_25/nz_256/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
         "CSDXMagFieldScanAr/nout_2_timestep_25/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
        )

    tSlices = (slice(80,240), slice(80,210))

    scanParameter = "B0"

    convertToPhysical = True

    xInd            = None
    yInd            = 16
    tSlice          = None
    nPoints         = 1
    equallySpace    = "x"
    steadyStatePath = steadyStatePaths[0]

    useMultiProcess = False

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

    print("\n\nTesting growth rates driver")
    dGR = DriverPhaseShift(
                     # DriverPhaseShift
                     dmp_folders                          ,\
                     scanCollectPaths                     ,\
                     steadyStatePaths                     ,\
                     tSlices                              ,\
                     scanParameter                        ,\
                     indicesArgs                          ,\
                     indicesKwargs                        ,\
                     plotSuperKwargs                      ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     useMultiProcess = useMultiProcess,\
                          )
    dGR.driverPhaseShift()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
