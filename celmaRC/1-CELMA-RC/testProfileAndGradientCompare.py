#!/usr/bin/env python

"""Post-processor test for profileAndGradientCompare"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.profAndGradCompare import driverProfAndGradCompare

#{{{profAndGradCompareTest
def profAndGradCompareTest():
    """
    Runs the comparison of profile and gradients test
    """

    collectPaths =\
        (\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/",\
         "CSDXMagFieldScanAr/nout_5000_timestep_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_1/"\
        )

    steadyStatePath = "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"

    varName = "n"

    savePath = "."
    convertToPhysical = True
    xSlice = 16
    ySlice = 16
    zSlice = 0
    tSlice = None

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting the comparison of profile and gradients")
    driverProfAndGradCompare(varName          ,\
                             collectPaths     ,\
                             steadyStatePath  ,\
                             convertToPhysical,\
                             xSlice           ,\
                             ySlice           ,\
                             zSlice           ,\
                             tSlice           ,\
                             plotSuperKwargs  ,\
                            )
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
