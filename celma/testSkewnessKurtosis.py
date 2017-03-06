#!/usr/bin/env python

"""Post-processor test for skewness and kurtosis"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.skewnessKurtosis import DriverSkewnessKurtosis, driverSkewnessKurtosis

#{{{skewKurtTest
def skewKurtTest():
    """
    Runs the test for the skewness and kurtosis
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

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting the skewness and kurtosis")
    driverSkewnessKurtosis(collectPaths     ,\
                           varName          ,\
                           convertToPhysical,\
                           mode             ,\
                           indicesArgs      ,\
                           indicesKwargs    ,\
                           plotSuperKwargs  ,\
                          )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the test for the skewness and kurtosis
    """

    dmp_folders  = ("CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1",)
    collectPaths =\
       (\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
       )

    useMultiProcess = False

    varName           = "n"
    convertToPhysical = True
    mode              = "fluct"

    yInd              = 16
    zInd              = 0
    tSlice            = None

    indicesArgs   = (None, yInd, zInd)
    indicesKwargs = {"tSlice" : tSlice}

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : None,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting the skewness and kurtosis driver")
    dSK = DriverSkewnessKurtosis(
                     # DriverSkewnessKurtosiss
                     dmp_folders                   ,\
                     indicesArgs                   ,\
                     indicesKwargs                 ,\
                     plotSuperKwargs               ,\
                     varName           = varName   ,\
                     mode              = mode      ,\
                     # DriverPointsSuperClass
                     convertToPhysical = convertToPhysical,\
                     # DriverSuperClass
                     collectPaths  = collectPaths ,\
                     useMultiProcess = useMultiProcess,\
                          )
    dSK.driverSkewnessKurtosis()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
