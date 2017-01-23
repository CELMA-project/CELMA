#!/usr/bin/env python

"""Post-processor test for blobs"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.blobs import DriverBlobs, prepareBlobs, driverWaitingTimePulse

#{{{testWaitingTimePulse
def testWaitingTimePulse():
    """
    Runs the test for the blobs
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

    print("\n\nTesting blobs")
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
    Runs the test for the blobs
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

    print("\n\nTesting blobs")
    ccb = prepareBlobs(collectPaths     ,\
                       slices           ,\
                       pctPadding       ,\
                       convertToPhysical,\
                      )
    driverBlobTimeTrace(ccb, plotSuperKwargs)
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
                        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/",\
       )

    xInd   = 26
    yInd   = 16
    zInd   = 0
    tSlice = slice(1200,None)
    slices = (xInd, yInd, zInd, tSlice)

    pctPadding = 400

    convertToPhysical = True

    useSubProcess = False

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
                     convertToPhysical,\
                     plotSuperKwargs  ,\
                     pctPadding       ,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useSubProcess = useSubProcess,\
                    )
    dB.driverBlobs()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    blobsTest()
#    driverTest()
