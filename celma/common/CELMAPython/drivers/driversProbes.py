#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import PerpPlaneProbes, PlotProbes
import numpy as np

#{{{DriversProbes
class DriversProbes(StatsAndSignalsDrivers):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self                          ,\
                 *args                         ,\
                 var                    = None ,\
                 yInd                   = None ,\
                 nProbes                = None ,\
                 steadyStatePath        = None ,\
                 maxMode                = 7    ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets additional member data.

        Parameters
        ----------
        *args : positional arguments
            See the constructor of StatsAndSignalsDrivers for details.
        var : str
            Variable to collect
        yInd : int
            yInd of the simulation data to use
        nProbes : int
            Number of probes to use
        steadyStatePath : str
            Path to the steady state
        maxMode : int
            Maximum modes to be plotted in FFT
        **kwargs : keyword arguments
            See the constructor of StatsAndSignalsDrivers and
            PlotProbes for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._var                    = var
        self._yInd                   = yInd
        self._nProbes                = nProbes
        self._maxMode                = maxMode

        if self._scanParameters:
            self._steadyStatePath =\
                    self._convertToCurrentScanParameters(self._steadyStatePath)
        else:
            self._steadyStatePath = steadyStatePath

        # Avgerage flux input
        self._uName     = "ExB"
        self._labelName = r"u_{E, \rho}"

        # Placeholder for the probes
        self._probes = None
    #}}}

    #{{{calcProbes
    def calcProbes(self):
        """ Calculates the statistics of the probes """
        # Create the probes
        self._probes = PerpPlaneProbes(\
                      self._var                                            ,\
                      paths                  = self._paths                 ,\
                      yInd                   = self._yInd                  ,\
                      nProbes                = self._nProbes               ,\
                      convertToPhysical      = self._convertToPhysical     ,\
                      steadyStatePath        = self._steadyStatePath       ,\
                      radialProbesIndices    = None                        ,\
                     )

        # Create the probe
        self._probes.initializeInputOutput(self._probes.radialProbesIndices,\
                                     [self._probes.yInd],\
                                     [0])
        self._probes.calcStatsMoments()
        self._probes.calcPDFs()
        self._probes.calcPSDs()
        self._probes.calcAvgFluxThroughVolumeElement(self._probes.radialExB,\
                                                     self._uName)
        self._probes.calcFFTs()

        # Set the position of the FFT plot to be the middle of the probes
        keyPos = int(np.ceil(len(self._probes.probesKeys)/2))
        self._positionKey = self._probes.probesKeys[keyPos]
    #}}}

    #{{{plotProbes
    def plotProbes(self):
        """ Plot the statistics of the probes """

        # Calculate the probes if not already done
        if self._probes == None:
            self.calcProbes()

        # Create the probesPlotter
        probesPlotter = PlotProbes(\
                self._probes,\
                showPlot  = self._showPlot         ,\
                savePlot  = self._savePlot         ,\
                extension = self._extension        ,\
                savePath  = self._savePath         ,\
                pltSize   = self._pltSize          ,\
                                  )

        probesPlotter.plotTimeTrace()
        probesPlotter.plotPDFs()
        probesPlotter.plotPSDs()
        probesPlotter.plotAvgFluxThrougVolumeElement(\
                                            uName     = self._uName,\
                                            labelName = self._labelName)

        probesPlotter.plotZFFT(self._positionKey, self._maxMode)
    #}}}
#}}}
