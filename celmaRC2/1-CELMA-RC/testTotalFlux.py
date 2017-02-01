#!/usr/bin/env python

"""Post-processor test for totalFlux"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.totalFlux import DriverTotalFlux, driverTotalFlux

#{{{totalFluxTest
def totalFluxTest():
    """
    Runs the test for the total flux
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    xInd              = 16
    yInd              = 16
    tSlice            = None
    mode              = "normal"
    convertToPhysical = True

    savePath          = "."

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting total flux")
    driverTotalFlux(collectPaths     ,\
                    xInd             ,\
                    yInd             ,\
                    tSlice           ,\
                    mode             ,\
                    convertToPhysical,\
                    plotSuperKwargs  ,\
                   )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the total flux
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)
    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    useSubProcess     = False
    xInd              = None
    yInd              = None
    tSlice            = None
    mode              = "normal"
    convertToPhysical = True

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting total flux driver")
    dTF = DriverTotalFlux(\
                     # DriverTotalFlux
                     dmp_folders                          ,\
                     plotSuperKwargs                      ,\
                     xInd              = xInd             ,\
                     yInd              = yInd             ,\
                     tSlice            = tSlice           ,\
                     mode              = mode             ,\
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dTF.driverTotalFlux()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
