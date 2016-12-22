#!/usr/bin/env python

"""Post-processor test for field2D"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fields2D import (driver2DFieldPerpSingle,\
                              driver2DFieldParSingle,\
                              driver2DFieldPolSingle,\
                              driver2DFieldPerpParSingle,\
                             )

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

    savePath          = "."
    varName           = "lnN"
    convertToPhysical = True
    xSlice            = None
    ySlice            = 16
    zSlice            = None
    tSlice            = None
    fluct             = False
    varyMaxMin        = False

    print("\n\nTesting perp 2D field")
    driver2DFieldPerpSingle(collectPaths     ,\
                            savePath         ,\
                            varName          ,\
                            convertToPhysical,\
                            xSlice           ,\
                            ySlice           ,\
                            zSlice           ,\
                            tSlice           ,\
                            fluct            ,\
                            varyMaxMin       ,\
                            xguards  = False ,\
                            yguards  = False ,\
                            showPlot = False ,\
                            savePlot = True  ,\
                                            )
    print("Success!\n\n")
#}}}

#{{{singleParTest
def singleParTest():
    """
    Runs the single par test for field2D
    """

    collectPaths =\
        (\
        "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_7.8633_geom_Ly_275.2144_input_B0_0.1_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/",\
        )

    savePath          = "."
    varName           = "jPar"
    convertToPhysical = True
    xSlice            = None
    ySlice            = None
    zSlice            = 0
    tSlice            = None
    fluct             = True
    varyMaxMin        = False

    print("\n\nTesting par 2D field")
    driver2DFieldParSingle(collectPaths     ,\
                            savePath         ,\
                            varName          ,\
                            convertToPhysical,\
                            xSlice           ,\
                            ySlice           ,\
                            zSlice           ,\
                            tSlice           ,\
                            fluct            ,\
                            varyMaxMin       ,\
                            xguards  = False ,\
                            yguards  = False ,\
                            showPlot = False ,\
                            savePlot = True  ,\
                                            )
    print("Success!\n\n")
#}}}

#{{{singlePolTest
def singlePolTest():
    """
    Runs the single pol test for field2D
    """

    collectPaths =\
        (\
        "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
        )

    savePath          = "."
    varName           = "uEPar"
    convertToPhysical = True
    xSlice            = 16
    ySlice            = None
    zSlice            = None
    tSlice            = None
    fluct             = False
    varyMaxMin        = False

    print("\n\nTesting pol 2D field")
    driver2DFieldPolSingle(collectPaths      ,\
                            savePath         ,\
                            varName          ,\
                            convertToPhysical,\
                            xSlice           ,\
                            ySlice           ,\
                            zSlice           ,\
                            tSlice           ,\
                            fluct            ,\
                            varyMaxMin       ,\
                            xguards  = False ,\
                            yguards  = False ,\
                            showPlot = False ,\
                            savePlot = True  ,\
                                             )
    print("Success!\n\n")
#}}}

#{{{singlePerpParTest
def singlePerpParTest():
    """
    Runs the single perpPar test for field2D
    """
    collectPaths =\
        (\
        "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/",\
       )

    savePath          = "."
    varName           = "uIPar"
    convertToPhysical = True
    xSlice            = None
    ySlice            = None
    zSlice            = None
    yInd              = 50
    zInd              = 1
    tSlice            = None
    fluct             = False
    varyMaxMin        = False

    print("\n\nTesting perpPar 2D field")
    driver2DFieldPerpParSingle(collectPaths     ,\
                               savePath         ,\
                               varName          ,\
                               convertToPhysical,\
                               xSlice           ,\
                               ySlice           ,\
                               zSlice           ,\
                               yInd             ,\
                               zInd             ,\
                               tSlice           ,\
                               fluct            ,\
                               varyMaxMin       ,\
                               xguards  = False ,\
                               yguards  = False ,\
                               showPlot = False ,\
                               savePlot = True  ,\
                                            )
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    singlePerpParTest()
