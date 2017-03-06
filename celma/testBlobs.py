#!/usr/bin/env python

"""Post-processor test for blobs"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.blobs import (DriverBlobs           ,\
                           driverBlobTimeTraces  ,\
                           driverPlot2DData      ,\
                           driverWaitingTimePulse,\
                           prepareBlobs          ,\
                          )

#{{{testWaitingTimePulse
def testWaitingTimePulse():
    """
    Runs the waiting time test for the blobs
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1000,None)
    slices = (xInd, yInd, zInd, tSlice)
    pctPadding = 400

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

    print("\n\nTesting waiting time")
    ccb = prepareBlobs(collectPaths     ,\
                       slices           ,\
                       pctPadding       ,\
                       convertToPhysical,\
                      )
    driverWaitingTimePulse(ccb, plotSuperKwargs)
    print("Success!\n\n")
#}}}

#{{{testBlobTimeTrace
def testBlobTimeTrace():
    """
    Runs the time trace test for the blobs
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1000,None)
    slices = (xInd, yInd, zInd, tSlice)
    pctPadding = 400

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

    print("\n\nTesting blob time trace")
    ccb = prepareBlobs(collectPaths     ,\
                       slices           ,\
                       pctPadding       ,\
                       convertToPhysical,\
                      )
    driverBlobTimeTraces(ccb, plotSuperKwargs)
    print("Success!\n\n")
#}}}

#{{{test2DDataPerp
def test2DDataPerp():
    """
    Runs the 2D perp test for the blobs
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1000,None)
    slices = (xInd, yInd, zInd, tSlice)
    pctPadding = 400

    mode  = "perp"
    fluct = True

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

    print("\n\nTesting blobs 2D perp")
    ccb = prepareBlobs(collectPaths     ,\
                       slices           ,\
                       pctPadding       ,\
                       convertToPhysical,\
                      )
    driverPlot2DData(ccb, mode, fluct, plotSuperKwargs)
    print("Success!\n\n")
#}}}

#{{{test2DDataPar
def test2DDataPar():
    """
    Runs the 2D par test for the blobs
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1000,None)
    slices = (xInd, yInd, zInd, tSlice)
    pctPadding = 400

    mode  = "par"
    fluct = False

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

    print("\n\nTesting blobs 2D par")
    ccb = prepareBlobs(collectPaths     ,\
                       slices           ,\
                       pctPadding       ,\
                       convertToPhysical,\
                      )
    driverPlot2DData(ccb, mode, fluct, plotSuperKwargs)
    print("Success!\n\n")
#}}}

#{{{test2DDataPol
def test2DDataPol():
    """
    Runs the 2D pol test for the blobs
    """

    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1000,None)
    slices = (xInd, yInd, zInd, tSlice)
    pctPadding = 400

    mode  = "pol"
    fluct = True

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

    print("\n\nTesting blobs 2D pol")
    ccb = prepareBlobs(collectPaths     ,\
                       slices           ,\
                       pctPadding       ,\
                       convertToPhysical,\
                      )
    driverPlot2DData(ccb, mode, fluct, plotSuperKwargs)
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the blobs
    """

    dmp_folders  = (    "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)

    collectPaths =\
       (\
                        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
                        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1200,None)
    slices = (xInd, yInd, zInd, tSlice)

    pctPadding = 400
    normed     = False
    convertToPhysical = True
    useMultiProcess = False

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }


    print("\n\nTesting the blobs driver")
    dB = DriverBlobs(\
                     # DriverBlobs
                     dmp_folders      ,\
                     slices           ,\
                     pctPadding       ,\
                     convertToPhysical,\
                     plotSuperKwargs  ,\
                     normed = normed  ,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                    )
    dB.driverAll()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
