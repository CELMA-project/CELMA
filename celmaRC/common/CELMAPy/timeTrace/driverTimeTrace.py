#!/usr/bin/env python

"""
Contains drivers for the time traces
"""

from ..commonDrivers import CommonPostProcessingDriver
from ..superClasses import PointsSuperClass
from .plotTimeTrace import PlotTimeTrace
from .calcTimeTrace import calcTimeTrace

#{{{DriverTimeTrace
class DriverTimeTrace(PointsSuperClass, CommonPostProcessingDriver):
    """
    Class which handles the time trace data.
    """

    #{{{Constructor
    def __init__(self,\
                 *args,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets the member data

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
        self._timeTrace = None
    #}}}

    #{{{getTimeTraces
    def getTimeTraces(self):
        """Obtain the timeTrace"""
        # Create the probes
        self._timeTrace = calcTimeTrace(\
                self._paths,\
                self._varName,\
                self._xInd,\
                self._yInd,\
                self._zInd,\
                convertToPhysical = self._convertToPhysical,\
                mode              = self._mode,\
                tSlice            = self._tSlice,\
                )
    #}}}

    #{{{plotTimeTrace
    def plotTimeTrace(self):
        """Plots the timeTrace"""

        # Calculate the probes if not already done
        if self._timeTrace is None:
            self.getTimeTraces()

        # Create the timeTracePlotter
        timeTracePlotter = PlotTimeTrace(\
                self._paths                                ,\
                self._timeTrace                            ,\
                convertToPhysical = self._convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        timeTracePlotter.plotTimeTrace()
    #}}}
#}}}
