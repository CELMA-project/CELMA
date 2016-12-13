#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import collectEnergy, PlotEnergy

#{{{DriverTimeTrace
class DriverTimeTrace(CommonPostProcessingDriver):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self,\
                 *args,\
                 paths,\
                 xInd,\
                 yInd,\
                 zInd,\
                 tSlice,\
                 varName   = "n",\
                 timeTrace = None,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets path and pltSize

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        paths  : tuple of strings
            The paths to collect from
        xInd : [int|sequence of ints]
            xInd to use when collecting.
            If given as a sequence, the length of the sequence must
            match the other input dimensions.
        yInd : [int|sequence of ints]
            yInd to use when collecting.
            If given as a sequence, the length of the sequence must
            match the other input dimensions.
        zInd : [int|sequence of ints]
            zInd to use when collecting.
            If given as a sequence, the length of the sequence must
            match the other input dimensions.
        tSlice : [slice|sequence of slices]
            Slice of t to use when collecting.
            If given as a sequence, the length of the sequence must
            match the other input dimensions.
        varName : str
            Name to collect.
        timeTrace : array
            Array of the time trace.
            Will be collected if not given.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._paths   = paths
        self._xInd    = xInd
        self._yInd    = yInd
        self._zInd    = zInd
        self._tSlice  = tSlice
        self._varName = varName

        self._pltName = getPltName(varName)
        self._pltSize = (12, 9)

        # Placeholder for the timeTrace
        self._timeTrace = None
    #}}}

    #{{{getTimeTraces
    def getTimeTraces(self):
        """Obtain the timeTrace"""
        # Create the probes
        self._timeTrace = calcTimeTrace(self._paths,\
                                        self._varName,\
                                        self._xInd,\
                                        self._yInd,\
                                        self._zInd,\
                                        self._tSlice,\
                                        self._convertToPhysical)
    #}}}

    #{{{plotTimeTrace
    def plotTimeTrace(self):
        """Plots the timeTrace"""

        # Calculate the probes if not already done
        if self._timeTrace == None:
            self.getTimeTraces()

        # Create the energyPlotter
        energyPlotter = PlotEnergy(\
                self._paths                                ,\
                self._energy                               ,\
                convertToPhysical = self._convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        energyPlotter.plotKinEnergy("ions")
        energyPlotter.plotKinEnergy("electrons")
        energyPlotter.plotPotEnergy()
    #}}}
#}}}
