#!/usr/bin/env python

"""
Contains drivers for the radial flux
"""

from ..commonDrivers import CommonPostProcessingDriver
from ..superClasses import PointsSuperClass
from .plotRadialFlux import PlotRadialFlux
from .calcRadialFlux import calcRadialFlux

#{{{DriverRadialFlux
class DriverRadialFlux(PointsSuperClass, CommonPostProcessingDriver):
    """
    Class which handles the radial flux data.
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
        self._RadialFlux = None
    #}}}

    #{{{getRadialFlux
    def getRadialFlux(self):
        """Obtain the RadialFlux"""
        # Create the probes
        self._RadialFlux = calcRadialFlux(self._paths,\
                            self._varName,\
                            self._xInd,\
                            self._yInd,\
                            self._zInd,\
                            converToPhysical = self.convertToPhysical,\
                            mode             = self._mode,\
                            tSlice           = self._tSlice,\
                            )
    #}}}

    #{{{plotRadialFlux
    def plotRadialFlux(self):
        """Plots the radial flux"""

        # Calculate the probes if not already done
        if self._timeTrace == None:
            self.getTimeTraces()

        # Create the energyPlotter
        RadialFluxPlotter = PlotRadialFlux(\
                self._paths                                ,\
                self._RadialFlux                           ,\
                convertToPhysical = self.convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        RadialFluxPlotter.plotRadialFlux()
    #}}}
#}}}
