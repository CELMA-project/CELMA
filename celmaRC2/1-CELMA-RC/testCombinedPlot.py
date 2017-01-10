#!/usr/bin/env python

"""Post-processor test for combinedPlots"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.combinedPlots import DriverCombinedPlots, driverCombinedPlots

#{{{combinedPlotTest
def combinedPlotTest():
    """
    Runs the test for the combined plots
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    steadyStatePath   = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    varName           = "n"
    convertToPhysical = True
    mode              = "fluct"

    yInd              = 16
    zInd              = 128
    tSlice            = None

    savePath          = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting combined plot")
    driverCombinedPlots(collectPaths     ,\
                        steadyStatePath  ,\
                        varName          ,\
                        convertToPhysical,\
                        mode             ,\
                        yInd             ,\
                        zInd             ,\
                        tSlice           ,\
                        plotSuperKwargs  ,\
                       )

    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the combined plots
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)
    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    steadyStatePath   = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    varName           = "n"
    convertToPhysical = True
    mode              = "fluct"

    yInd              = 16
    zInd              = 128
    tSlice            = slice(10, None, None)

    savePath          = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }


    useSubProcess = False

    print("\n\nTesting combined plot driver")
    dTT = DriverCombinedPlots(
                     # DriverCombinedPlots
                     dmp_folders                ,\
                     steadyStatePath            ,\
                     yInd                       ,\
                     zInd                       ,\
                     tSlice                     ,\
                     plotSuperKwargs            ,\
                     varName           = varName,\
                     mode              = mode   ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dTT.driverCombinedPlots()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
