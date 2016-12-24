#!/usr/bin/env python

"""
Contains single driver and driver class for the time traces
"""

from ..superClasses import CollectAndCalcPointsSuperClass
from .collectAndCalcTimeTrace import collectAndCalcTimeTrace
from .plotTimeTrace import PlotTimeTrace
import os

#{{{driverTimeTraceSingle
def driverTimeTraceSingle(collectPaths     ,\
                          fieldPlotType    ,\
                          convertToPhysical,\
                          xInd             ,\
                          yInd             ,\
                          zInd             ,\
                          tSlice           ,\
                          equallySpace     ,\
                          steadyStatePath  ,\
                          mode             ,\
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
    fieldPlotType : str
        What predefined fieldPlotType to plot.
        See getCollectFieldsAndPlotOrder for details.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    xInd : [None|int|sequence of ints]
        xInd to use when collecting.
            * If None: This constructor will use the index of the
                       largest gradient in "n" from the steady state
                       path. This value will be the center-index of
                       nPoints with equidistant spacing
            * If int: If yInd and zInd are given as sequences:
                      The same as None, but the center index will be
                      given by the input int.
                      Else: The value will be copied nPoints times
            * If sequence: All the indices are given
        In all cases, the resulting length of the tuple must match
        the other physical dimensions.
    yInd : [int|sequence of ints]
        The same as xInd (except the None possibility), but for the y-index.
    zInd : [int|sequence of ints]
        The same as xInd (except the None possibility), but for the z-index.
    tSlice : [None|sequence of slices]
        If given this is the  slice of t to use when collecting.
        The length of the sequence must match the other input
        dimensions.
    nPoints : [None|int]
        Size of the sequence. Ignored if xInd, yInd and zInd are all
        sequences.
        See xInd for details.
    equallySpace : ["x", "y", "z"]
        If there is any ambiguity of which coordinate to equally
        space around one value, this variable will be used.
        Default is "x".
    steadyStatePath: str
        Path to find the gradient in "n". Only used if xInd is None
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class
    """
    #}}}


    # Create collect object
    cctt = CollectAndCalcTimeTrace(collectPaths                         ,\
                                   mode              = mode             ,\
                                   convertToPhysical = convertToPhysical,\
                                   )

    # Set the slice
    cctt.setIndices(xInd, yInd, zInd, tSlice,\
                    nPoints, equallySpac, steadyStatePath)

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

YOU ARE HERE
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
