#!/usr/bin/env python

"""Post-processor test for profileAndGradientCompare"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.radialProfile import (DriverRadialProfile,\
                                   driverProfAndGradCompare,\
                                   driverPosOfFluct)

#{{{profAndGradCompareTest
def profAndGradCompareTest():
    """
    Runs the comparison of profile and gradients test
    """

    collectPaths =\
        (\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/",\
        )

    steadyStatePath = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    varName = "n"

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

    print("\n\nTesting the comparison of profile and gradients")
    driverProfAndGradCompare(varName          ,\
                             collectPaths     ,\
                             steadyStatePath  ,\
                             convertToPhysical,\
                             yInd             ,\
                             plotSuperKwargs  ,\
                            )
    print("Success!\n\n")
#}}}

#{{{posOfFluctTest
def posOfFluctTest():
    """
    Runs the position of fluctation test
    """

    collectPaths =\
        (\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
        )

    steadyStatePath = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    var1Name = "n"
    var2Name = "phi"

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

    print("\n\nTesting the position of fluctuations")
    driverPosOfFluct(var1Name         ,\
                     var2Name         ,\
                     collectPaths     ,\
                     steadyStatePath  ,\
                     convertToPhysical,\
                     yInd             ,\
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

    useSubProcess = False

    varName           = "n"
    var2Name          = "phi"
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
    dRP = DriverRadialProfile(
                     # DriverRadialProfile
                     dmp_folders                          ,\
                     steadyStatePath                      ,\
                     yInd                                 ,\
                     tSlice                               ,\
                     plotSuperKwargs                      ,\
                     varName           = varName          ,\
                     var2Name          = var2Name         ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dRP.driverPosOfFluct()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
