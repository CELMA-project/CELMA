#!/usr/bin/env python

"""Post-processor test for magnitudeSpectrum"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.magnitudeSpectrum import DriverMagnitudeSpectrum, driverMagnitudeSpectrum

#{{{magnitudeSpectrumTest
def magnitudeSpectrumTest():
    """
    Runs the test for the magnitude spectrum
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    varName           = "n"
    convertToPhysical = True

    xInd            = 16
    yInd            = 16
    tSlice          = slice(200,400)
    nPoints         = 1
    equallySpace    = "x"
    steadyStatePath = None

    indicesArgs   = (xInd, yInd)
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

    print("\n\nTesting magnitude spectrum")
    driverMagnitudeSpectrum(collectPaths     ,\
                            varName          ,\
                            convertToPhysical,\
                            indicesArgs      ,\
                            indicesKwargs    ,\
                            plotSuperKwargs  ,\
                           )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the magnitude spectrum
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)
    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    useSubProcess = False

    varName           = "n"
    convertToPhysical = True
    nModes            = 7

    xInd              = None
    yInd              = 16
    zInd              = None
    tSlice            = slice(1000,None)
    nPoints           = 3
    equallySpace      = "x"

    steadyStatePath   = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

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

    print("\n\nTesting magnitude spectrum driver")
    dMS = DriverMagnitudeSpectrum(
                     # DriverMagnitudeSpectrum
                     dmp_folders              ,\
                     indicesArgs              ,\
                     indicesKwargs            ,\
                     plotSuperKwargs          ,\
                     varName         = varName,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                          )
    dMS.driverMagnitudeSpectrum()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
