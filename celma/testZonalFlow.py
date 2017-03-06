#!/usr/bin/env python

"""Post-processor test for zonalFlow"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.zonalFlow import (DriverZonalFlow,\
                               driverZonalFlow)

#{{{zonalFlowTest
def zonalFlowTest():
    """
    Runs the zonal flow test
    """

    collectPaths =\
        (\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
        )

    steadyStatePath = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    savePath = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    convertToPhysical = True
    yInd = 16
    tSlice = None

    print("\n\nTesting the position of fluctuations")
    driverZonalFlow(collectPaths     ,\
                    steadyStatePath  ,\
                    convertToPhysical,\
                    yInd             ,\
                    tSlice           ,\
                    plotSuperKwargs  ,\
                   )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the time traces
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)

    collectPaths =\
        (\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
        )

    steadyStatePath = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    useMultiProcess     = False
    convertToPhysical = True
    yInd              = 16
    tSlice            = None

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting radial profile driver")
    dRP = DriverZonalFlow(
                     # DriverZonalFlow
                     dmp_folders                          ,\
                     steadyStatePath                      ,\
                     yInd                                 ,\
                     tSlice                               ,\
                     plotSuperKwargs                      ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dRP.driverZonalFlow()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
