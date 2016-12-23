#!/usr/bin/env python

"""
Contains drivers for the time traces
"""

from ..superClasses import CollectAndCalcPointsSuperClass
from .collectAndCalcTimeTrace import collectAndCalcTimeTrace
from .plotTimeTrace import PlotTimeTrace
import os

#{{{DriverTimeTrace
class DriverTimeTrace(CollectAndCalcPointsSuperClass):
    """
    Class which handles the time trace data.
    """

    #{{{Constructor
    def __init__(self                  ,\
                 *args                 ,\
                 timeStampFolder = True,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets the member data
        * Updates the savePath and makes the folder

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        timeStampFolder : bool
            Whether or not to timestamp the folder
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

        # Update the savePath
        firstPathPart = os.path.dirname(self._savePath)
        if timeStampFolder:
            timePath = os.path.basename(self._savePath)
        else:
            timePath = ""
        self._savePath = os.path.join(firstPathPart, "timeTraces", timePath)

        # Make dir if not exists
        if not os.path.exists(self._savePath):
            os.makedirs(self._savePath)
    #}}}

    #{{{getTimeTraces
    def getTimeTraces(self):
        """Obtain the timeTrace"""
        # Create the probes
        self._timeTrace = collectAndCalcTimeTrace(\
                self._collectPaths                         ,\
                self._varName                              ,\
                self._xInd                                 ,\
                self._yInd                                 ,\
                self._zInd                                 ,\
                convertToPhysical = self.convertToPhysical,\
                mode              = self._mode             ,\
                tSlice            = self._tSlice           ,\
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
                convertToPhysical = self.convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        timeTracePlotter.plotTimeTrace()
    #}}}
#}}}
