#!/usr/bin/env python
"""
Contains single driver and driver class for the phase shifts
"""
from ..superClasses import DriverPointsSuperClass
from .collectAndCalcAnalyticGrowthRates import CollectAndCalcAnalyticGrowthRates
from .plotGrowthRates import PlotGrowthRates
from multiprocessing import Process

#{{{driverPhaseShift
def driverPhaseShift(steadyStatePaths,\
                              scanParameter   ,\
                              yInd            ,\
                              plotSuperKwargs ,\
                             ):
    #{{{docstring
    """
    Driver for plotting phase shifts.

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

    import pdb; pdb.set_trace()
    print(paramDataFrame)

    # Plot
    pagr = PlotGrowthRates(uc         ,\
                           **plotSuperKwargs)
    pagr.setData("", analyticalGRDataFrame, positionTuple, analytic=True)
    pagr.plotSaveShowGrowthRates()
#}}}

#{{{DriverPhaseShift
class DriverPhaseShift(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the phase shifts.
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

    #{{{driverPhaseShift
    def driverPhaseShift(self):
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
        if self._useSubProcess:
            processes = Process(target = driverPhaseShift, args = args)
            processes.start()
        else:
            driverPhaseShift(*args)
    #}}}
#}}}
