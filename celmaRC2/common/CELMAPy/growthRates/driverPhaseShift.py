#!/usr/bin/env python
"""
Contains single driver and driver class for the phase shifts
"""
from ..superClasses import DriverPointsSuperClass
from .collectAndCalcPhaseShift import CollectAndCalcPhaseShift
from .driverGrowthRates import DriverGrowthRates
from .plotPhaseShift import PlotPhaseShift
from multiprocessing import Process

#{{{driverPhaseShift
def driverPhaseShift(collectArgs    ,\
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

    ccps = CollectAndCalcPhaseShift()
    ccps.setCollectArgs(*collectArgs)
    ccps.setGetDataArgs(*getDataArgs[1:-1])
    apsdf, psdf, positionTuple, uc = ccps.getData()

    # Plot
    paps = PlotPhaseShift(uc         ,\
                          **plotSuperKwargs)
    paps.setData(apsdf, psdf, positionTuple)
    paps.plotSaveShowPhaseShift()
#}}}

#{{{DriverPhaseShift
class DriverPhaseShift(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the phase shifts.
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
        self._collectArgs=DriverGrowthRates.makeCollectArgs(scanCollectPaths,\
                                                            steadyStatePaths,\
                                                            tSlices         ,\
                                                            scanParameter)
        self._getDataArgs=DriverGrowthRates.makeGetDataArgs(""              ,\
                                                           convertToPhysical,\
                                                           indicesArgs      ,\
                                                           indicesKwargs    ,\
                                                           0)

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"phaseShift"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverPhaseShift
    def driverPhaseShift(self):
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
        if self._useMultiProcess:
            processes = Process(target = driverPhaseShift, args = args)
            processes.start()
        else:
            driverPhaseShift(*args)
    #}}}
#}}}
