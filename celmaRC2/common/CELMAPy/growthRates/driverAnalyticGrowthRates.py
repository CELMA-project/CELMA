#!/usr/bin/env python
"""
Contains single driver and driver class for the analytic growth rates
"""
from ..superClasses import DriverPointsSuperClass
from .collectAndCalcAnalyticGrowthRates import CollectAndCalcAnalyticGrowthRates
from .plotGrowthRates import PlotGrowthRates
from multiprocessing import Process

#{{{driverAnalyticGrowthRates
def driverAnalyticGrowthRates(steadyStatePaths,\
                              scanParameter   ,\
                              yInd            ,\
                              plotSuperKwargs ,\
                             ):
    #{{{docstring
    """
    Driver for plotting analytic growth rates.

    Parameters
    ----------
    steadyStatePaths : tuple
        Paths to the steady states.
    scanParameter : str
        String of the parameter which is scanned.
    yInd : int
        The y-index to use.
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

    print(paramDataFrame)

    # Plot
    pagr = PlotGrowthRates(uc         ,\
                           **plotSuperKwargs)
    pagr.setData("", analyticalGRDataFrame, positionTuple, analytic=True)
    pagr.plotSaveShowGrowthRates()
#}}}

#{{{DriverAnalyticGrowthRates
class DriverAnalyticGrowthRates(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the analytic growth rates.
    """

    #{{{Constructor
    def __init__(self            ,\
                 dmp_folders     ,\
                 steadyStatePaths,\
                 scanParameter   ,\
                 yInd            ,\
                 plotSuperKwargs ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Set the member data
            * Updates the plotSuperKwargs

        Parameters
        ----------
        dmp_folders : tuple
            Tuple of the dmp_folder (output from bout_runners).
        steadyStatePaths : tuple
            Paths to the steady states.
        scanParameter : str
            String of the parameter which is scanned.
        yInd : int
            The y-index to use.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._steadyStatePaths = steadyStatePaths
        self._scanParameter    = scanParameter
        self._yInd             = yInd
        self._plotSuperKwargs  = plotSuperKwargs

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"growthRates"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverAnalyticGrowthRates
    def driverAnalyticGrowthRates(self):
        #{{{docstring
        """
        Wrapper to driverFourierMode
        """
        #}}}
        args =  (\
                self._steadyStatePaths,\
                self._scanParameter   ,\
                self._yInd            ,\
                self._plotSuperKwargs ,\
               )
        if self._useMultiProcess:
            processes = Process(target = driverAnalyticGrowthRates, args = args)
            processes.start()
        else:
            driverAnalyticGrowthRates(*args)
    #}}}
#}}}
