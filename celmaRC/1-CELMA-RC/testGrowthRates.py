#!/usr/bin/env python

"""Post-processor test for growthRates"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.growthRates import DriverGrowthRates, driverGrowthRates

#{{{growthRatesTest
def growthRatesTest():
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
         "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
         "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
        )

    startInds = (\
                 80,\
                 80,\
                )
    endInds   = (\
                 240,\
                 210,\
                )

    scanParameter = "B0"

    collectArgs = DriverGrowthRates.makeCollectArgs(scanCollectPaths,\
                                                    steadyStatePaths,\
                                                    startInds       ,\
                                                    endInds         ,\
                                                    scanParameter)

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd            = 16
    yInd            = 16
    tSlice          = None
    nPoints         = 3
    equallySpace    = "x"
    steadyStatePath = None

    indicesArgs   = (xInd, yInd)
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

    print("\n\nTesting growth rates")
    driverGrowthRates(collectArgs    ,\
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
    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    useSubProcess = False

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd              = None
    yInd              = 16
    zInd              = None
    tSlice            = None
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

    print("\n\nTesting growth rates driver")
    dFM = DriverGrowthRates(
                     # DriverGrowthRates
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

if __name__ == "__main__":
    growthRatesTest()
