#!/usr/bin/env python

"""Post-processor test for field1D"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fields1D import Driver1DFields, driver1DFieldSingle

#{{{singleParallelTest
def singleParallelTest():
    """
    Runs the single parallel test for field1D
    """
    collectPaths =\
        (\
          "CSDXMagFieldScanAr/nout_2_timestep_2000.0/nz_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/",\
          "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"\
        )

    savePath = "."
    fieldPlotType = "mainFields"
    convertToPhysical = True
    xSlice = 0
    ySlice = None
    zSlice = 0
    tSlice = None
    mode   = "parallel"
    hyperIncluded = False

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }


    print("\n\nTesting parallel 1D field")
    driver1DFieldSingle(collectPaths     ,\
                        fieldPlotType    ,\
                        convertToPhysical,\
                        xSlice           ,\
                        ySlice           ,\
                        zSlice           ,\
                        tSlice           ,\
                        mode             ,\
                        hyperIncluded    ,\
                        plotSuperKwargs  ,\
                        )
    print("Success!\n\n")
#}}}

#{{{singleRadialTest
def singleRadialTest():
    """
    Runs the single radial test for field1D
    """
    collectPaths =\
        (\
          "CSDXMagFieldScanAr/nout_2_timestep_2000.0/nz_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/",\
          "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"\
        )

    savePath = "."
    fieldPlotType = "mainFields"
    convertToPhysical = True
    xSlice = None
    ySlice = 16
    zSlice = 0
    tSlice = None
    mode   = "radial"
    hyperIncluded = False

    plotSuperKwargs = {\
                        "showPlot"     : False,\
                        "savePlot"     : True,\
                        "savePath"     : savePath,\
                        "savePathFunc" : None,\
                        "extension"    : None,\
                        "dmp_folders"  : None,\
                       }

    print("\n\nTesting perpendicular 1D")
    driver1DFieldSingle(collectPaths     ,\
                        fieldPlotType    ,\
                        convertToPhysical,\
                        xSlice           ,\
                        ySlice           ,\
                        zSlice           ,\
                        tSlice           ,\
                        mode             ,\
                        hyperIncluded    ,\
                        plotSuperKwargs  ,\
                        )
    print("Success!\n\n")
#}}}

#{{{driverTest
def driverTest():
    """
    Runs the single test for field1D
    """
    dmp_folders  = ("CSDXMagFieldScanAr/nout_2_timestep_2000.0/nz_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0",)
    collectPaths =\
        (\
          "CSDXMagFieldScanAr/nout_2_timestep_2000.0/nz_1/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/",\
          "CSDXMagFieldScanAr/nout_2_timestep_50/nz_256/geom_Lx_4.718_geom_Ly_165.1286_input_B0_0.06_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-1-expand_0/"\
        )

    useMultiProcess = False

    convertToPhysical = True
    hyperIncluded = False

    xSlice = None
    ySlice = None
    zSlice = None
    tSlice = None
    xInd = 0
    yInd = 16
    zInd = 0
    guardSlicesAndIndicesKwargs = {\
                                   "xguards" : False ,\
                                   "yguards" : False ,\
                                   "xSlice"  : xSlice,\
                                   "ySlice"  : ySlice,\
                                   "zSlice"  : zSlice,\
                                   "tSlice"  : tSlice,\
                                   "xInd"    : xInd  ,\
                                   "yInd"    : yInd  ,\
                                   "zInd"    : zInd  ,\
                                  }
    # Plot super kwargs
    theRunName = "test"
    savePathFunc = "scanWTagSaveFunc"
    plotSuperKwargs = {\
                        "showPlot"       : False       ,\
                        "savePlot"       : True        ,\
                        "savePath"       : None        ,\
                        "savePathFunc"   : savePathFunc,\
                        "extension"      : None        ,\
                        "dmp_folders"    : None        ,\
                        "timeStampFolder": True        ,\
                        "theRunName"     : theRunName  ,\
                       }

    print("\n\nTesting 1D driver")
    d1DF = Driver1DFields(\
                   # Driver1DFields
                   dmp_folders                    ,\
                   plotSuperKwargs                ,\
                   guardSlicesAndIndicesKwargs    ,\
                   boussinesq      = False        ,\
                   hyperIncluded   = hyperIncluded,\
                   # DriverPlotFieldsSuperClass
                   convertToPhysical = convertToPhysical,\
                   # DriverSuperClass
                   collectPaths  = collectPaths ,\
                   useMultiProcess = useMultiProcess,\
                  )

    d1DF.driver1DFieldsAll()
    print("Success!\n\n")
#}}}

if __name__ == "__main__":
    driverTest()
