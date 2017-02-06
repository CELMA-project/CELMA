#!/usr/bin/env python

"""Post-processor test for the analytic growthRates"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.growthRates import DriverAnalyticGrowthRates, driverAnalyticGrowthRates

#{{{analyticGrowthRatesTest
def analyticGrowthRatesTest():
    """
    Runs the test for the analytic growth rates
    """

    steadyStatePaths =\
        (\
         "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
         "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
        )

    scanParameter = "B0"
    yInd          = 16

    savePath          = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting analytic growth rates")
    driverAnalyticGrowthRates(steadyStatePaths,\
                              scanParameter   ,\
                              yInd            ,\
                              plotSuperKwargs ,\
                             )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the analytic growth rates
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)

    steadyStatePaths =\
        (\
         "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_1.5727_geom_Ly_55.0429_input_B0_0.02_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
         "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
        )

    scanParameter = "B0"
    yInd          = 16
    useMultiProcess = False

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting analytic growth rates driver")
    dAGR = DriverAnalyticGrowthRates(
                                     # DriverAnalyticGrowthRates
                                     dmp_folders     ,\
                                     steadyStatePaths,\
                                     scanParameter   ,\
                                     yInd            ,\
                                     plotSuperKwargs ,\
                                     # DriverSuperClass
                                     useMultiProcess = useMultiProcess,\
                                    )
    dAGR.driverAnalyticGrowthRates()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
