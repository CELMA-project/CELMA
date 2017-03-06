#!/usr/bin/env python

"""Post-processor for growthRates"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.growthRates import DriverAnalyticGrowthRates

#{{{analyticGrowthRatesPlot
def analyticGrowthRatesPlot(dmp_folders     ,\
                            steadyStatePaths,\
                            scanParameter   ,\
                            plotSuperKwargs ,\
                           ):
    #{{{docstring
    """
    Runs the standard analytic growth rates plot

    Parameters
    ----------
    dmp_folders : tuple
        Tuple of the dmp_folder (output from bout_runners).
    steadyStatePaths : tuple
        Paths to the steady states.
    scanParameter : str
        String of the parameter which is scanned.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    **kwargs : keyword arguments
        See parent class for details.
    """
    #}}}

    yInd          = 16
    useMultiProcess = False

    dAGR = DriverAnalyticGrowthRates(
                                     # DriverAnalyticGrowthRates
                                     dmp_folders     ,\
                                     steadyStatePaths,\
                                     scanParameter   ,\
                                     yInd            ,\
                                     plotSuperKwargs ,\
                                     # DriverSuperClass
                                     useMultiProcess = useMultiProcess,\
                                    )
    dAGR.driverAnalyticGrowthRates()
#}}}
