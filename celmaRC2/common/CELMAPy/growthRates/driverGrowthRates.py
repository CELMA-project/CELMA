#!/usr/bin/env python

"""
Contains single driver and driver class for the growth rates
"""

from ..superClasses import DriverPointsSuperClass
from .collectAndCalcGrowthRates import CollectAndCalcGrowthRates
from .plotGrowthRates import PlotGrowthRates
from multiprocessing import Process

#{{{driverGrowthRates
def driverGrowthRates(collectArgs    ,\
                      getDataArgs    ,\
                      plotSuperKwargs,\
                     ):
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
    ccgr = CollectAndCalcGrowthRates(*collectArgs)

    # Obtain the data
    growthRateDataFrame, positionTuple, uc =\
        ccgr.getData(*getDataArgs)

    # Plot
    pgr = PlotGrowthRates(uc         ,\
                          **plotSuperKwargs)
    pgr.setData(getDataArgs[0], growthRateDataFrame, positionTuple)
    pgr.plotSaveShowGrowthRates()
#}}}

#{{{DriverGrowthRates
class DriverGrowthRates(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the growth rates.
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

    @staticmethod
    #{{{makeCollectArgs
    def makeCollectArgs(scanCollectPaths,\
                        steadyStatePaths,\
                        tSlices         ,\
                        scanParameter):
        #{{{docstring
        """
        Sets the arguments to be used in the constructor of
        CollectAndCalcGrowthRates

        Parameters
        ----------
        See the constructor of DriverGrowthRates for details.

        Returns
        -------
        collectArgs : tuple
            Tuple containing arguments used in the collect function.
            See the constructor of CollectAndCalcGrowthRates for details.
        """
        #}}}

        collectArgs = (scanCollectPaths,\
                       steadyStatePaths,\
                       tSlices         ,\
                       scanParameter)

        return collectArgs
    #}}}

    @staticmethod
    #{{{makeGetDataArgs
    def makeGetDataArgs(varName          ,\
                        convertToPhysical,\
                        indicesArgs      ,\
                        indicesKwargs    ,\
                        nModes):
        #{{{docstring
        """
        Sets the arguments to be used in
        CollectAndCalcGrowthRates.getData

        Parameters
        ----------
        See the constructor of DriverGrowthRates for details.

        Returns
        -------
        getDataArgs : tuple
            Tuple containing arguments used in
            CollectAndCalcGrowthRates.getData (including the keyword
            argument nModes).
            See the constructor of CollectAndCalcGrowthRates.getData for
            details.
        """
        #}}}

        getDataArgs = (varName          ,\
                       convertToPhysical,\
                       indicesArgs      ,\
                       indicesKwargs    ,\
                       nModes)

        return getDataArgs
    #}}}

    #{{{driverGrowthRates
    def driverGrowthRates(self):
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
