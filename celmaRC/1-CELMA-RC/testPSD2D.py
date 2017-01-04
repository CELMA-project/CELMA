#!/usr/bin/env python

"""Post-processor test for PSDs"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.PSD import DriverPSD, driverPSD2D

#{{{PSD2DTest
def PSD2DTest():
    """
    Runs the test for the 2D power spectral density
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    varName           = "n"
    convertToPhysical = False
    mode              = "fluct"

    yInd              = 16
    zInd              = 128
    tSlice            = None

    indicesArgs   = (None, yInd, zInd)
    indicesKwargs = {"tSlice" : tSlice}

    savePath = "."

    plotLimits = {"xlim":None,\
                  "ylim":None,\
                  "zlim":None}

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting 2D probability density function")
    driverPSD2D(collectPaths     ,\
                varName          ,\
                convertToPhysical,\
                mode             ,\
                indicesArgs      ,\
                indicesKwargs    ,\
                plotLimits       ,\
                plotSuperKwargs  ,\
               )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the 2D power spectral density
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)
    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    useSubProcess = False

    varName           = "n"
    convertToPhysical = True
    mode              = "fluct"

    yInd              = 16
    zInd              = 0
    tSlice            = None

    indicesArgs   = (None, yInd, zInd)
    indicesKwargs = {"tSlice" : tSlice}

    plotLimits = {"xlim":None      ,\
                  "ylim":(100, 3e4),\
                  "zlim":(-7,0)}

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting probability density function driver")
    dPSD = DriverPSD(
                     # DriverPSDs
                     dmp_folders                   ,\
                     indicesArgs                   ,\
                     indicesKwargs                 ,\
                     plotSuperKwargs               ,\
                     varName           = varName   ,\
                     mode              = mode      ,\
                     plotLimits        = plotLimits,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dPSD.driverPSD2D()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
