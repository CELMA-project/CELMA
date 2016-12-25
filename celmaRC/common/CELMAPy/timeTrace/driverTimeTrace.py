#!/usr/bin/env python

"""
Contains single driver and driver class for the time traces
"""

from ..superClasses import CollectAndCalcPointsSuperClass
from .collectAndCalcTimeTrace import CollectAndCalcTimeTrace
from .plotTimeTrace import PlotTimeTrace
import os

#{{{driverTimeTrace
def driverTimeTrace(collectPaths     ,\
                    varName          ,\
                    convertToPhysical,\
                    mode             ,\
                    indicesArgs      ,\
                    indicesKwargs    ,\
                    plotSuperKwargs  ,\
                    ):
    #{{{docstring
    """
    Driver for plotting time traces.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    indicesArgs : tuple
        Contains xInd, yInd and zInd.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    indicesArgs : dict
        Contains tslice, nPoints, equallySpace and steadyStatePath.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}


    # Create collect object
    cctt = CollectAndCalcTimeTrace(collectPaths                         ,\
                                   mode              = mode             ,\
                                   convertToPhysical = convertToPhysical,\
                                   )

    # Set the slice
    cctt.setIndices(*indicesArgs, **indicesKwargs)

    # Set name
    cctt.setVarName(varName)

    # Execute the collection
    tt = cctt.executeCollectAndCalc()
    tt = cctt.convertTo1D(tt)

    # Plot
    ptt = PlotTimeTrace(cctt.uc         ,\
                        **plotSuperKwargs)
    ptt.setData(tt, mode)
    ptt.plotSaveShowTimeTrace()
#}}}

# FIXME: YOU ARE HERE
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
        self._timeTrace = CollectAndCalcTimeTrace(\
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
