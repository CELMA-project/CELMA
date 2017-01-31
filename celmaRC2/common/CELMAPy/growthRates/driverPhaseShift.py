#!/usr/bin/env python
"""
Contains single driver and driver class for the phase shifts
"""
from ..superClasses import DriverPointsSuperClass
from .collectAndCalcAnalyticGrowthRates import CollectAndCalcAnalyticGrowthRates
from .collectAndCalcGrowthRates import CollectAndCalcGrowthRates
from .plotGrowthRates import PlotGrowthRates
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

# Collect the analytic phase shift
    # Create collect object
    # FIXME: Mixing of interfaces gives hard to read code
    steadyStatePaths = collectArgs[-3]
    scanParameter    = collectArgs[-1]
    yInd             = getDataArgs[2][1]
    ccagr = CollectAndCalcAnalyticGrowthRates(steadyStatePaths,\
                                              scanParameter,\
                                              yInd)

    # Obtain the data
    analyticalGRDataFrame, paramDataFrame, positionTuple, uc =\
        ccagr.getData()

# Do this for all the b0 scans
# Make a dict for the indices to extract from
    from scipy import signal
    # NOTE: The triangular window corresponds to the periodogram
    # estimate of the spectral density
    csd = signal.csd(n, phi, window="triang")

    # Find the magnitude
    # Expand the time dimension (only one point)
    csd = np.expand_dims(csd,axis=0)
    # Put into dict so in can be used in calcMagnitude
    # FIXME: Hack: Recast csd to a 2d signal, expand dimension so that
    # only one time, then can use
csdDictWMagnitues = CollectAndCalcGrowthRates.calcMagnitude(csdDict)
# Extract Magnitudes

    # Recast to 1d
    csd = csd[0,:]
find max index
avgPhaseShiftNPhi = np.angle(csd[maxInd])


    foo you are here
    import pdb; pdb.set_trace()
    print(paramDataFrame)

    # Plot
    pagr = PlotGrowthRates(uc         ,\
                           **plotSuperKwargs)
    pagr.setData("", analyticalGRDataFrame, positionTuple, analytic=True)
    pagr.plotSaveShowGrowthRates()
#}}}

# FIXME: NEW INPUTS
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
