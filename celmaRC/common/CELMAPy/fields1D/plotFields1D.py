#!/usr/bin/env python

"""
Contains functions for animating the 1D fields
"""

from ..superClasses import PlotAnim1DSuperClass
from ..plotHelpers import plotNumberFormatter
import os

#{{{PlotAnim1DRadial
class PlotAnim1DRadial(PlotAnim1DSuperClass):
    """
    Class for 1D radial plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (20,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot1DRadial

        * Calls the parent class
        * Sets the spatial title
        * Sets the xlabel

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        pltSize : tuple
            The size of the plot
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._pltSize = pltSize

        # Set the spatial part of the title
        self._spatTitle = "{}$,$ {}$,$"
        # Set the x-axis label
        self._xlabel = self._ph.rhoTxtDict["rhoTxtLabel"]
    #}}}

    #{{{setData
    def setData(self, radialDict, figName, savePath, plotOrder=None):
        #{{{docstring
        """
        Sets the radial data and set up the plotting

        Specifically this function will:
            * Set the data
            * Create the figures and axes (through a function call)
            * Set the colors (through a function call)
            * Set the plot order
            * Update the spatial title

        Parameters
        ----------
        radialDict : dict
            Dictionary of the data containing the following keys
                * var        - The collected variables.
                               A 2d array for each variable
                * "X"        - The abscissa of the variable
                * "time"     - The time trace
                * "zPos"     - The z position
                * "thetaPos" - The theta position
        figName : str
            Name of the figure
        savePath : str
            Destination to save the plot in.
        plotOrder : [None|sequence of str]
            If given: A sequence of the variable names in the order to
            plot them
        """
        #}}}

        self._X        = radialDict.pop("X")
        self._time     = radialDict.pop("time")
        self._vars     = radialDict

        zPos           = radialDict.pop("zPos")
        thetaPos       = radialDict.pop("thetaPos")

        self._figName  = figName
        self._savePath = savePath

        # Make axes and colors
        self._createFiguresAndAxes()
        self._setColors()

        # Set the plot order
        if plotOrder:
            self._plotOrder = tuple(plotOrder)
        else:
            self._plotOrder = tuple(sorted(list(self._vars.keys())))

        # Check for "ddt" in the variables
        self._ddtPresent = False

        for var in self._vars:
            if "ddt" in var:
                self._ddtVar = var
                self._ddtPresent = True
                # Remove ddt from the plotOrder
                self._plotOrder = list(self._plotOrder)
                self._plotOrder.remove(var)
                self._plotOrder = tuple(self._plotOrder)
                break

        # Update the title
        self._ph.zTxtDict["value"] = plotNumberFormatter(zPos, None)
        zTxt =\
            self._ph.zTxtDict["constZTxt"].format(self._ph.zTxtDict)
        self._ph.thetaTxtDict["value"] = plotNumberFormatter(thetaPos, None)
        thetaTxt =\
            self._ph.thetaTxtDict["constThetaTxt"].format(self._ph.thetaTxtDict)
        self._spatTitle = self._spatTitle.format(zTxt, thetaTxt)
    #}}}

    #{{{plotAndSaveProfile
    def plotAndSaveProfile(self):
        #{{{docstring
        """
        Performs the actual plotting of the radial plane
        """
        #}}}

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._figName, "radial", "1D"))

        # Initial plot
        self._initialRadialPlot()

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updateRadialAxInTime,\
                          len(self._time))
    #}}}

    #{{{_initialRadialPlot
    def _initialRadialPlot(self):
        #{{{docstring
        """
        Initial radial plot.

        The initial radial plot:
            * Calls the generic initial plot routine
            * Sets the title

        Subsequent animation is done with _updateRadialAxInTime.
        """
        #}}}

        # Initial plot the axes
        self._initialPlot()

        # Set the title
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[0], None)
        curTimeTxt = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}{}".format(self._spatTitle, curTimeTxt))
    #}}}

    #{{{_updateRadialAxInTime
    def _updateRadialAxInTime(self, tInd):
        #{{{docstring
        """
        Function which updates the data.

        Parameters
        ----------
        tInd : int
            The current t index.
        """
        #}}}

        for line, key in zip(self._lines, self._plotOrder):
            line.set_data(self._X, self._vars[key][tInd,:])

        # If ddt is present
        if self._ddtPresent:
            # NOTE: ddtLines is one longer than plotOrder
            for line, key in zip(self._ddtLines, self._plotOrder):
                line.set_data(self._X, self._vars[key][tInd,:])

            self._ddtLines[-1].set_data(self._X, self._vars[self._ddtVar][tInd,:])

        # Update the title
        self._ph.tTxtDict["value"] = plotNumberFormatter(self._time[tInd], None)
        curTimeTxt = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}{}".format(self._spatTitle, curTimeTxt))
    #}}}
#}}}

#{{{PlotAnim1DParallel
class PlotAnim1DParallel(PlotAnim1DSuperClass):
    """
    Class for 1D parallel plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (20,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot1DParallel

        * Calls the parent class
        * Sets the spatial title
        * Sets the xlabel

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        pltSize : tuple
            The size of the plot
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._pltSize = pltSize

        # Set the spatial part of the title
        self._spatTitle = "{}$,$ {}$,$"
        # Set the x-axis label
        self._xlabel = self._ph.rhoTxtDict["rhoTxtLabel"]
    #}}}

    #{{{setData
    def setData(self, parallelDict, figName, savePath, plotOrder=None):
        #{{{docstring
        """
        Sets the parallel data and set up the plotting

        Specifically this function will:
            * Set the data
            * Create the figures and axes (through a function call)
            * Set the colors (through a function call)
            * Set the plot order
            * Update the spatial title

        Parameters
        ----------
        parallelDict : dict
            Dictionary of the data containing the following keys
                * var        - The collected variables.
                               A 2d array for each variable
                * "X"        - The abscissa of the variable
                * "time"     - The time trace
                * "rhoPos"   - The rho position
                * "thetaPos" - The theta position
        figName : str
            Name of the figure
        savePath : str
            Destination to save the plot in.
        plotOrder : [None|sequence of str]
            If given: A sequence of the variable names in the order to
            plot them
        """
        #}}}

        self._X        = parallelDict.pop("X")
        self._time     = parallelDict.pop("time")
        self._vars     = parallelDict

        rhoPos         = parallelDict.pop("rhoPos")
        thetaPos       = parallelDict.pop("thetaPos")

        self._figName  = figName
        self._savePath = savePath

        # Make axes and colors
        self._createFiguresAndAxes()
        self._setColors()

        # Set the plot order
        if plotOrder:
            self._plotOrder = tuple(plotOrder)
        else:
            self._plotOrder = tuple(sorted(list(self._vars.keys())))

        # Check for "ddt" in the variables
        self._ddtPresent = False

        for var in self._vars:
            if "ddt" in var:
                self._ddtVar = var
                self._ddtPresent = True
                # Remove ddt from the plotOrder
                self._plotOrder = list(self._plotOrder)
                self._plotOrder.remove(var)
                self._plotOrder = tuple(self._plotOrder)
                break

        # Update the title
        self._ph.rhoTxtDict["value"] = plotNumberFormatter(rhoPos, None)
        rhoTxt =\
            self._ph.rhoTxtDict["constRhoTxt"].format(self._ph.rhoTxtDict)
        self._ph.thetaTxtDict["value"] = plotNumberFormatter(thetaPos, None)
        thetaTxt =\
            self._ph.thetaTxtDict["constThetaTxt"].format(self._ph.thetaTxtDict)
        self._spatTitle = self._spatTitle.format(rhoTxt, thetaTxt)
    #}}}

    #{{{plotAndSaveProfile
    def plotAndSaveProfile(self):
        #{{{docstring
        """
        Performs the actual plotting of the parallel plane
        """
        #}}}

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._figName, "parallel", "1D"))

        # Initial plot
        self._initialParallelPlot()

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updateParallelAxInTime,\
                          len(self._time))
    #}}}

    #{{{_initialParallelPlot
    def _initialParallelPlot(self):
        #{{{docstring
        """
        Initial parallel plot.

        The initial parallel plot:
            * Calls the generic initial plot routine
            * Sets the title

        Subsequent animation is done with _updateParallelAxInTime.
        """
        #}}}

        # Initial plot the axes
        self._initialPlot()

        # Set the title
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[0], None)
        curTimeTxt = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}{}".format(self._spatTitle, curTimeTxt))
    #}}}

    #{{{_updateParallelAxInTime
    def _updateParallelAxInTime(self, tInd):
        #{{{docstring
        """
        Function which updates the data.

        Parameters
        ----------
        tInd : int
            The current t index.
        """
        #}}}

        for line, key in zip(self._lines, self._plotOrder):
            line.set_data(self._X, self._vars[key][tInd,:])

        # If ddt is present
        if self._ddtPresent:
            # NOTE: ddtLines is one longer than plotOrder
            for line, key in zip(self._ddtLines, self._plotOrder):
                line.set_data(self._X, self._vars[key][tInd,:])

            self._ddtLines[-1].set_data(self._X, self._vars[self._ddtVar][tInd,:])

        # Update the title
        self._ph.tTxtDict["value"] = plotNumberFormatter(self._time[tInd], None)
        curTimeTxt = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}{}".format(self._spatTitle, curTimeTxt))
    #}}}
#}}}
