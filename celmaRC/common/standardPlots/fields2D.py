#!/usr/bin/env python

"""Post-processor for fields2D animation"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.fields2D import Driver2DFields

#{{{fields2DAnimation
def fields2DAnimation(dmp_folders, collectPaths, steadyStatePath,\
                      varName = "n", fluct=False):
    #{{{docstring
    """
    Runs the standard fields2D animation

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folders.
    collectPaths : tuple
        Tuple of the paths to collect from.
    steadyStatePath : str
        Path to the steady state.
    varName : str
        Variable to animate.
    fluct : bool
        Whether or not to plot the fluctuations.
    """
    #}}}

    useSubProcess = False

    convertToPhysical = True
    tSlice            = None
    varyMaxMin = True if fluct else False

    # DriverPlotFieldsSuperClass
    xSlice = None
    ySlice = None
    zSlice = None
    tSlice = None
    xInd   = None
    yInd   = 16
    zInd   = 0
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
    plotSuperKwargs = {\
                        "showPlot"        : False      ,\
                        "savePlot"        : True       ,\
                        "savePath"        : None       ,\
                        "savePathFunc"    : None       ,\
                        "extension"       : None       ,\
                        "dmp_folders"     : dmp_folders,\
                        "timeStampFolder" : False      ,\
                       }

    d2DF = Driver2DFields(
                   # Driver2DFields
                   dmp_folders                   ,\
                   plotSuperKwargs               ,\
                   guardSlicesAndIndicesKwargs   ,\
                   varName         = varName     ,\
                   fluct           = fluct       ,\
                   varyMaxMin      = varyMaxMin  ,\
                   # DriverPlotFieldsSuperClass
                   convertToPhysical = convertToPhysical,\
                   # DriverSuperClass
                   collectPaths  = collectPaths ,\
                   useSubProcess = useSubProcess,\
                  )

    if fluct:
        d2DF.setXIndToMaxGradInN(steadyStatePath)
        d2DF.driver2DFieldsPerpPol()
    else:
        d2DF.driver2DFieldsPerpPar()
#}}}
