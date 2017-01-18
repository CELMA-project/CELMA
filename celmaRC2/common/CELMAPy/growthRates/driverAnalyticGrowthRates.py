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
        Paths to the steady states
    scanParameter : str
        String of the parameter which is scanned
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

    # Plot
    pagr = PlotGrowthRates(uc         ,\
                           **plotSuperKwargs)
    pagr.setData("", analyticalGRDataFrame, positionTuple, analytic=True)
    pagr.plotSaveShowGrowthRates()
#}}}

# FIXME:
#{{{DriverAnalyticGrowthRates
class DriverAnalyticGrowthRates(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the analytic growth rates.
    """

    #{{{Constructor
    def __init__(self                    ,\
                 dmp_folders             ,\
                 scanCollectPaths        ,\
                 steadyStatePaths        ,\
                 tSlices                 ,\
                 scanParameter           ,\
                 indicesArgs             ,\
                 indicesKwargs           ,\
                 plotSuperKwargs         ,\
                 varName           = "n" ,\
                 nModes            = 7   ,\
                 convertToPhysical = True,\
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
        scanCollectPaths : tuple of tuple of strings
            One tuple of strings for each value of the scan which will
            be used in collective collect.
        steadyStatePaths : tuple
            Path to the steady state simulation.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        tSlices : tuple of slices
            The time slices to use for each folder in the scan order.
            This must manually be selected to the linear phase found in
            from the fourier moedes.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        scanParameter : str
            String segment representing the scan
        indicesArgs : tuple
            Tuple of indices to use when collecting.
            NOTE: Only one spatial point should be used.
        indicesKwargs : dict
            Keyword arguments to use when setting the indices for
            collecting.
            NOTE: Only one spatial point should be used.
        varName : str
            Name of variable to find the growth rates of.
        nModes : int
            Number of modes.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._collectArgs = self.makeCollectArgs(scanCollectPaths,\
                                                 steadyStatePaths,\
                                                 tSlices         ,\
                                                 scanParameter)
        self._getDataArgs = self.makeGetDataArgs(varName          ,\
                                                 convertToPhysical,\
                                                 indicesArgs      ,\
                                                 indicesKwargs    ,\
                                                 nModes)

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
                 self._collectArgs    ,\
                 self._getDataArgs    ,\
                 self._plotSuperKwargs,\
                )
        if self._useSubProcess:
            processes = Process(target = driverGrowthRates, args = args)
            processes.start()
        else:
            driverGrowthRates(*args)
    #}}}
#}}}
