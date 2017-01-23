#!/usr/bin/env python

"""Class for blobs plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotBlobTimeTrace
class PlotBlobTimeTrace(PlotSuperClass):
    """
    Class which contains the time traces together with the plotting
    configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (15,10), **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets member data

        Parameters
        ----------
        pltSize : tuple
            The size of the plot
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, theDict, blobOrHole):
        #{{{docstring
        """
        Sets the waiting time and pulse widths to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        theAverage : dict
            Dictionary of the averaged blob/hole.
            Contains the keys "n", "time" and "pos".
        blobOrHole : ["blobs"|"holes"]
            Whether blob data or hole data is given as an input.
        """
        #}}}

        # Set the member data
        self._var        = theDict.pop("n")
        self._time       = theDict.pop("time")
        self._pos        = theDict.pop("pos")
        self._blobOrHole = blobOrHole

        self._avg = True if self._time[0] < 0 else False

        avg = "-avg" if self._avg else ""

        fileName = "timeTrace-{}{}".format(self._blobOrHole, avg)

        self._prepareLabels()

        # Set the fileName
        fileName = os.path.join(self._savePath,\
            "{}-{}".format(, self._blobOrHole))

        # Get the number
        files = glob(fileName + "*")
        if len(files) !=0:
            files = sorted(files)
            nr    = int(files[-1].split("-")[-1].split(".")[0])+1
        else:
            nr = 0

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}-{}.{}".format(fileName, nr, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Set x and y labels
        if self.uc.convertToPhysical:
            units  = r"$[\mathrm{s}]$"
            yUnits = r"$[\mathrm{m}^{-3}]$"
            normalization  = ""
            yNormalization = ""
        else:
            units  = ""
            yUnits = ""
            normalization  = r"$\om_{ci}$"
            yNormalization = r"$/n_0$"

        self._xLabel = r"$t$ " + normalization + units
        self._yLabel = r"$\widetilde{n}$ " +\
                       yNormalization +\
                       yUnits

        # Set title
        if self._avg:
            self._title = r"$\mathrm{{Average \;{}}}$".format(self._blobOrHole)
        else:
            self._title = r"$\mathrm{{{}}}$".\
                    format(self._blobOrHole.capitalize())
    #}}}

    #{{{plotSaveShowTimeTrace
    def plotSaveShowTimeTrace(self):
        """
        Performs the actual plotting.
        """

        # Create the plot
        fig, ax = plt.subplots(figsize = self._pltSize)

        # Waiting time
        ax.plot(self._var ,\
                self._time,\
               )
        ax.set_xlabel(self._xLabel)
        ax.set_ylabel(self._yLabelW)
        ax.set_title(self._wTitle)

        # Set the title
        fig.suptitle(self._title)

        # Make the plot look nice
        self._ph.makePlotPretty(ax, rotation = 45, legend = False)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}

#{{{PlotTemporalStats
class PlotTemporalStats(PlotSuperClass):
    """
    Class which contains the waiting time and pulse width data together
    with the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (30,10), bins = "sqrt", **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets member data

        Parameters
        ----------
        pltSize : tuple
            The size of the plot
        bins : [int|str]
            The binning, see numpy.histogram for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._pltSize = pltSize
        self._bins    = bins
    #}}}

    #{{{setData
    def setData(self, waitingTimes, pulseWidths, blobOrHole, normed=False):
        #{{{docstring
        """
        Sets the waiting time and pulse widths to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        waitingTimes : tuple
            Tuple containing the waiting times.
        pulseWidths : tuple
            Tuple containing the waiting times.
        blobOrHole : ["blobs"|"holes"]
            Whether blob data or hole data is given as an input.
        normed : bool
            Wheter or not the histogram should be normed
        """
        #}}}

        # Set the member data
        self._waitingTimes = waitingTimes
        self._pulseWidths  = pulseWidths
        self._blobOrHole   = blobOrHole
        self._normed       = normed

        self._prepareLabels()

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
            "{}-{}".format("temporalStats", self._blobOrHole))

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Set x and y labels
        if self.uc.convertToPhysical:
            units  = r"$[\mathrm{s}]$"
            yUnits = r"$[\mathrm{s}^{-1}]$"
            normalization  = ""
            yNormalization = ""
        else:
            units  = ""
            yUnits = ""
            normalization  = r"$\om_{ci}$"
            yNormalization = r"$/\om_{ci}$"

        self._xLabel  = r"$t$ " + normalization + units
        if self._normed:
            self._yLabelW = r"$\mathrm{PDF}(\mathrm{Waiting\;time}$ " +\
                            yNormalization +\
                            "$) $" +\
                            yUnits
            self._yLabelP = r"$\mathrm{PDF}(\mathrm{Pulse\;width}$ " +\
                            yNormalization +\
                            "$) $" +\
                            yUnits
        else:
            self._yLabelW = r"$\mathrm{Waiting\;time\;count}$"
            self._yLabelP = r"$\mathrm{Pulse\;width\;count}$"

        # Set title
        self._title = r"$\mathrm{{Temporal \;statistics \;for \;{}}}$".\
                        format(self._blobOrHole)
        self._wTitle = r"$\mathrm{Waiting \;time}$"
        self._pTitle = r"$\mathrm{Pulse \;width}$"
    #}}}

    #{{{plotSaveShowTemporalStats
    def plotSaveShowTemporalStats(self):
        """
        Performs the actual plotting.
        """

        # Create the plot
        fig, (wAx, pAx) =\
                plt.subplots(nrows=1, ncols=2, figsize = self._pltSize)

        # Waiting time
        wAx.hist(self._pulseWidths    ,\
                 self._bins           ,\
                 normed = self._normed,\
                 alpha  = 0.75)
        wAx.set_xlabel(self._xLabel)
        wAx.set_ylabel(self._yLabelW)
        wAx.set_title(self._wTitle)

        # Pulse width
        pAx.hist(self._waitingTimes   ,\
                 self._bins           ,\
                 normed = self._normed,\
                 alpha  = 0.75)
        pAx.set_xlabel(self._xLabel)
        pAx.set_ylabel(self._yLabelP)
        pAx.set_title(self._pTitle)

        # Set the title
        fig.suptitle(self._title)

        # Make the plot look nice
        self._ph.makePlotPretty(wAx, rotation = 45, legend = False)
        self._ph.makePlotPretty(pAx, rotation = 45, legend = False)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.2, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
