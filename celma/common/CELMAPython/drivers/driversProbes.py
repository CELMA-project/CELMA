#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import PerpPlaneProbes, PlotProbes
from multiprocessing import Process
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
                 tIndSaturatedTurb      = None ,\
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
        self._var               = var
        self._tIndSaturatedTurb = tIndSaturatedTurb
        self._yInd              = yInd
        self._nProbes           = nProbes
        self._maxMode           = maxMode

        if self._scanParameters:
            self._steadyStatePath =\
                    self._convertToCurrentScanParameters(steadyStatePath)
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
                      self._var                                    ,\
                      paths               = self._paths            ,\
                      yInd                = self._yInd             ,\
                      nProbes             = self._nProbes          ,\
                      convertToPhysical   = self._convertToPhysical,\
                      steadyStatePath     = self._steadyStatePath  ,\
                      tIndSaturatedTurb   = self._tIndSaturatedTurb,\
                      radialProbesIndices = None                   ,\
                     )

        # Create the probe
        self._probes.initializeInputOutput(self._probes.radialProbesIndices,\
                                     [self._probes.yInd],\
                                     [0])

        # FIXME: Use subprocesses to calculate
        #        May have to define a calling classe, see
        #        http://stackoverflow.com/questions/35717109/python-class-object-sharing-between-processes-created-using-multiprocessing-modu
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
