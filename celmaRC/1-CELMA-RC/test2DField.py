#!/usr/bin/env python

"""Post-processor test for field2D"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fields2D import driver2DFieldSingle

#{{{singlePerpTest
def singlePerpTest():
    """
    Runs the single perp test for field2D
    """
    collectPaths =\
        (\
        "CSDXMagFieldScanAr/nout_1000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-2-linearPhase1_0/",\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
       )

    savePath = "."
    varName = "lnN"
    convertToPhysical = True
    xSlice = None
    ySlice = 16
    zSlice = None
    tSlice = None
    fluct  = False
    mode   = "perp"

    print("\n\nTesting perp 2D field")
    driver2DFieldSingle(collectPaths     ,\
                        savePath         ,\
                        varName          ,\
                        convertToPhysical,\
                        xSlice           ,\
                        ySlice           ,\
                        zSlice           ,\
                        tSlice           ,\
                        fluct            ,\
                        mode             ,\
                        xguards  = False ,\
                        yguards  = False ,\
                        showPlot = False ,\
                        savePlot = True  ,\
                                        )
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    singlePerpTest()
