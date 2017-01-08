#!/usr/bin/env python

"""Post-processor for fields1D animation"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fields1D import Driver1DFields

#{{{fields1DAnimation
def fields1DAnimation(dmp_folders, collectPaths, hyperIncluded=False):
    #{{{docstring
    """
    Runs the standard fields1D animation

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders
    collectPaths : tuple
        Tuple of the paths to collect from
    hyperIncluded : bool
        If hyper viscosities are used.
    """
    #}}}

    useSubProcess     = True
    convertToPhysical = True

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

    plotSuperKwargs = {\
                        "showPlot"        : False ,\
                        "savePlot"        : True  ,\
                        "savePath"        : None  ,\
                        "savePathFunc"    : None  ,\
                        "extension"       : None  ,\
                        "dmp_folders"     : None  ,\
                        "timeStampFolder" : False ,\
                       }

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
                   useSubProcess = useSubProcess,\
                  )

    d1DF.driver1DFieldsAll()
#}}}
