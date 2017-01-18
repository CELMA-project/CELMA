#!/usr/bin/env python
"""
Contains single driver and driver class for the analytic growth rates
"""

from .collectAndCalcAnalyticGrowthRates import CollectAndCalcAnalyticGrowthRates
# from multiprocessing import Process

#{{{driverGrowthRates
def driverGrowthRates(steadyStatePaths,\
                      scanParameter   ,\
                      yInd            ,\
                      # plotSuperKwargs,\
                     ):
# FIXME
    #{{{docstring
    """
    Driver for plotting growth rates.

    Parameters
    ----------
    collectArgs : tuple
        Tuple containing arguments used in the collect function.
        See the constructor of CollectAndCalcGrowthRates for details.
    getDataArgs : tuple
        Tuple containing arguments used in
        CollectAndCalcGrowthRates.getData (including the keyword
        argument nModes).
        See the constructor of CollectAndCalcGrowthRates.getData for
        details.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    # Create collect object
    ccagr = CollectAndCalcAnalyticGrowthRates(steadyStatePaths,\
                                              scanParameter,\
                                              yInd)

    # Obtain the data
    analyticalGRDataFrame, paramDataFrame, positionTuple, uc =\
        ccagr.getData()

    import pdb; pdb.set_trace()
    a = 1
#}}}
