#!/usr/bin/env python

"""
Contains drivers for calculation of the growth rates
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import PerpPlaneProbes, calcGrowthRate, plotGrowthRates
from multiprocessing import Process

#{{{DriversGrowthRates
class DriversGrowthRates(StatsAndSignalsDrivers):
    """
    Class which handles calculation and plotting of the growth rates
    """

    #{{{Constructor
    def __init__(self                   ,\
                 *args                  ,\
                 var              = None,\
                 scanParam        = None,\
                 yInd             = None,\
                 steadyStatePaths = None,\
                 maxMode          = 7   ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets additional member data.

        NOTE: The order of steadyStatePaths must agree with the order of
              self._paths

        Parameters
        ----------
        *args : positional arguments
            See the constructor of StatsAndSignalsDrivers for details.
        var : str
            Variable to collect
        scanParam : str
            Parameter which is scanned
        yInd : int
            yInd of the simulation data to use
        steadyStatePaths : iterable
            Path to the steady states.
            NOTE! The order of steadyStatePaths must agree with the
                  order of self._paths.
        maxMode : int
            Maximum modes to make a dispersion relation diagram from
        **kwargs : keyword arguments
            See the constructor of StatsAndSignalsDrivers and
FIXME:
            PlotGrowthRates for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._var              = var
        self._scanParam        = scanParam
        self._yInd             = yInd
        self._maxMode          = maxMode
        self._steadyStatePaths = steadyStatePaths

        # Placeholder for the growthRates
        self._growthRatesDict = None
    #}}}

    #{{{calcProbes
    def calcProbes(self):
        """ Calculates the statistics of the probes """

        # Make the growthRatesDict a dictionary
        self._growthRatesDict = {}

        # Create the probes (one for each scan point)
        # Sort the folders (ensures that they come in the right order)
        self._dmp_folder.sort()
        self._steadyStatePaths.sort()
        for path, steadyStatePath in zip(self._paths, self._steadyStatePaths):
            # FIXME: The loop is parallelizable
            # We only find the growth rates at position of highest gradient,
            # as we from theory expect the highest growth rates to be there
            curProbe = PerpPlaneProbes(\
                          self._var                                    ,\
                          paths               = path                   ,\
                          yInd                = self._yInd             ,\
                          convertToPhysical   = self._convertToPhysical,\
                          steadyStatePath     = steadyStatePath        ,\
                          nProbes             = 1                      ,\
                          radialProbesIndices = None                   ,\
                         )

            # Initialize the probes in order to calculate
            curProbe.initializeInputOutput(curProbe.radialProbesIndices,\
                                           [curProbe.yInd],\
                                           [0])

            # Calculate the FFT
            curProbe.calcFFTs()

            # Get the positionKey of the probe
            positionKey = list(curProbe.results.keys())[0]

            # Obtain the current scan value
            if hasattr(path, "__iter__") and type(path) != str:
                # The path is a list
                scanValue = path[0].split("_")
            else:
                # The path is a string
                scanValue = path.split("_")
            # +1 as the value is immediately after the scan parameter
            scanValue = scanValue[scanValue.index(self._scanParam)+1]
            # The time and the time is clipped with three to clip away
            # where initial modes decay
            initClip = 3
            # We will also clip so that approximately only the linear
            # mode is present (less data to process)
            linClip = curProbe.results[positionKey]["zFFTLinearIndex"]
            if linClip <= initClip:
                message = ("{0}{1}WARNING: "\
                           "Could not find a proper startpoint for the "
                           "linear stage{1}{0}")
                linClip = None
                print(message.format("\n"*2, "!"*5))
            # Clipping modes and time
            modes = curProbe.results[positionKey]["zFFT"][initClip:linClip, :]
            time  = curProbe.time[initClip:linClip]

            # Find the growth rates
            self._growthRatesDict["{}={}".format(self._scanParam, scanValue)]=\
                    calcGrowthRate(modes   = modes,\
                                   time    = time ,\
                                   maxMode = self._maxMode)
    #}}}

    #{{{plotGrowthRates
    def plotGrowthRates(self):
        """ Plot the growth rates"""

        # Calculate the probes if not already done
        if self._growthRatesDict == None:
            self.calcProbes()

# FIXME: Make one plot per mode
        # Create the probesPlotter
        probesPlotter = plotGrowthRates()
        #PlotProbes(\
        #        self._probes,\
        #        showPlot  = self._showPlot         ,\
        #        savePlot  = self._savePlot         ,\
        #        extension = self._extension        ,\
        #        savePath  = self._savePath         ,\
        #        pltSize   = self._pltSize          ,\
        #                          )

        if self._useSubProcess:
            #{{{ Function call through subprocess
            proc = {}
            # Create process
            proc["plotTimeTrace"] =\
                    Process(\
                            target = probesPlotter.plotTimeTrace,\
                            args   = ()                         ,\
                            kwargs = {}
                           )
            proc["plotPDFs"] =\
                    Process(\
                            target = probesPlotter.plotPDFs,\
                            args   = ()                    ,\
                            kwargs = {}
                           )
            proc["plotPSDs"] =\
                    Process(\
                            target = probesPlotter.plotPSDs,\
                            args   = ()                    ,\
                            kwargs = {}
                           )
            proc["plotAvgFluxThroughVolumeElement"] =\
                    Process(\
                target =  probesPlotter.plotAvgFluxThroughVolumeElement,\
                args   = (self._uName, self._labelName)                ,\
                kwargs = {}
                           )
            proc["plotZFFT"] =\
                    Process(\
                            target = probesPlotter.plotZFFT            ,\
                            args   = (self._positionKey, self._maxMode),\
                            kwargs = {}
                           )

            # Start processes
            for key in proc.keys():
                proc[key].start()
            # Wait for process to finish
            for key in proc.keys():
                proc[key].join()
            #}}}
        else:
            #{{{Normal function call
                probesPlotter.plotTimeTrace()
                probesPlotter.plotPDFs()
                probesPlotter.plotPSDs()
                probesPlotter.plotAvgFluxThroughVolumeElement(\
                                                    uName     = self._uName,\
                                                    labelName = self._labelName)
                probesPlotter.plotZFFT(self._positionKey, self._maxMode)
            #}}}
    #}}}
#}}}
