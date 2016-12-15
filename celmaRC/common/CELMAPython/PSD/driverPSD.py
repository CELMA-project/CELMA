#!/usr/bin/env python

"""
Contains drivers for the PSD
"""

from ..commonDrivers import CommonPostProcessingDriver
from ..superClasses import PointsSuperClass
from .plotPSD import PlotPSD
from .calcPSD import calcPSD

#{{{DriverPSD
class DriverPSD(PointsSuperClass, CommonPostProcessingDriver):
    """
    Class which handles the PSD data.
    """

    #{{{Constructor
    def __init__(self,\
                 *args,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets pltSize

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._pltSize = (12, 9)

        # Placeholder for the timeTrace
        self._PSD = None
    #}}}

    #{{{getPSD
    def getPSD(self):
        """Obtain the PSD"""
        # Create the probes
        self._PSD = calcPSD(self._paths,\
                            self._varName,\
                            self._xInd,\
                            self._yInd,\
                            self._zInd,\
                            converToPhysical = self._convertToPhysical,\
                            mode             = self._mode,\
                            tSlice           = self._tSlice,\
                            )
    #}}}

    #{{{plotPSD
    def plotPSD(self):
        """Plots the PSD"""

        # Calculate the probes if not already done
        if self._timeTrace == None:
            self.getTimeTraces()

        # Create the energyPlotter
        PSDPlotter = PlotPSD(\
                self._paths                                ,\
                self._PSD                                  ,\
                convertToPhysical = self._convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        PSDPlotter.plotPSD()
    #}}}
#}}}
