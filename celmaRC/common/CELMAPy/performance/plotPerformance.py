#!/usr/bin/env python

"""Class for performance plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, qualCMap
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotPerformance
class PlotPerformance(PlotSuperClass):
    """
    Class which contains the performance data and the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (15,10), **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor

        Parameters
        ----------
        pltSize : tuple
            The size of the plot
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the plot size
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, performance, allFolders=False):
        #{{{docstring
        """
        Sets the performance to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        performance : dict
            Although the keys may vary depending on the settings, the
            standard keys are:
                * time      - The normalized time in the simulation
                * RHSPrTime - Number of right hand side evaluations before a
                              time step per time ste:
                * WallTime  - Physical time used on the step
                * Calc      - Percentage of time used on arithmetical
                              calculation
                * Inv       - Percentage of time used on laplace inversion
                * Comm      - Percentage of time used on communication
                * I/O       - Percentage of time used on inupt/output
                * SOLVER    - Percentage of time used in the solver
            The values are stored in numpy arrays.
        allFolders : bool
            If true, the fileName will be tagged with allFolders.
        """
        #}}}

        # Set the member data
        self._performance = performance
        keys = performance.keys()

        self._colors = qualCMap(np.linspace(0, 1, len(keys)))

        self._prepareLabels()

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath, "performance")

        # Whether or not all folders has been used
        if allFolders:
            self._fileName += "allFolders"

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Prepare legends
        self._legends = []
        self._legends.append("$\mathrm{Arithmetics}$")
        self._legends.append("$\mathrm{Laplace \quad inversions}$")
        self._legends.append("$\mathrm{Communication}$")
        self._legends.append("$\mathrm{Input/output}$")
        self._legends.append("$\mathrm{Time \quad solver}$")
        self._legends = sorted(self._legends)
        self._legends = tuple(self._legends)

        # Set the y-labels
        self._RHSEvalYlabel = r"$\mathrm{RHS \quad evaluations}/\mathrm{Time step}$"
        self._pctYlabel     = r"$\mathrm{Percentage \quad} [\%]$"

        # Set the time label
        self._timeLabel = self._ph.tTxtDict["tTxtLabel"]
    #}}}

    #{{{plotSaveShowPerformance
    def plotSaveShowPerformance(self):
        """
        Performs the actual plotting.
        """

        # Create the plot
        fig, (RHSAx, pctAx) =\
                plt.subplots(nrows=2, figsize=self._pltSize, sharex=True)

        # Pop RHS evals and time as treated separately
        RHSEvals = self._performance.pop("RHSPrTime")
        time     = self._performance.pop("time")

        # RHS evals axis
        RHSAx.plot(time, RHSEvals, color=self._colors[0])
        RHSAx.set_ylabel(self._RHSEvalYlabel)

        # Percentages axis
        keys = sorted(self._performance.keys())
        for key, color, label in zip(keys, self._colors[1:], self._legends):
            pctAx.plot(time,\
                       self._performance[key],\
                       color=color, label=label, alpha=0.7)

        # Set axis labels
        pctAx.set_ylabel(self._pctYlabel)
        pctAx.set_xlabel(self._timeLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(RHSAx, rotation = 45, legend=False)
        self._ph.makePlotPretty(pctAx, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
