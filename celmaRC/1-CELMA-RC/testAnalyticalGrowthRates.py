#!/usr/bin/env python

"""Post-processor test for the analytic growthRates"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

# FIXME
# from CELMAPy.growthRates import DriverGrowthRates, driverGrowthRates
from CELMAPy.growthRates import driverGrowthRates

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

    print("\n\nTesting analytic growth rates")
    driverGrowthRates(steadyStatePaths,\
                          scanParameter   ,\
                          yInd            ,\
                          # plotSuperKwargs,\
                     )
    print("Success!\n\n")
#}}}