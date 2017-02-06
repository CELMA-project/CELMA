#!/usr/bin/env python

"""Class for blobs plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import SizeMaker
from glob import glob
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import os

#{{{PlotBlobOrHoleTimeTraceSingle
class PlotBlobOrHoleTimeTraceSingle(PlotSuperClass):
    """
    Class which contains the time traces together with the plotting
    configuration of either blobs, or holes.
    """

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)
    #}}}

    #{{{setData
    def setData(self, theDict, blobOrHole):
        #{{{docstring
        """
        Sets the waiting time and pulse widths to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        theDict : dict
            Dictionary of the blob/hole.
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

        self._title = self._prepareLabels(self._blobOrHole)

        self._setFileName(fileName)
    #}}}

    #{{{_setFileName
    def _setFileName(self, fileName):
        #{{{docstring
        """
        Sets self._fileName
        """
        #}}}
        # Set the fileName
        fileName = os.path.join(self._savePath, fileName)

        # Get the number
        files = glob(fileName + "*")
        if len(files) != 0:
            nrs = sorted([int(curFile.split("-")[-1].split(".")[0])\
                          for curFile in files])
            nr  = nrs[-1] + 1
        else:
            nr = 0

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}-{}.{}".format(fileName, nr, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self, blobOrHole):
        #{{{docstring
        """
        Prepares the labels for plotting.

        Parameters
        ----------
        blobOrHole : ["blob"|"hole"]
            Whether blobs or holes are being considered.

        Returns
        -------
        axTitle : str
            The axis title
        """
        #}}}

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

        self._xLabel = r"$t$ " if self._time[0] >= 0 else r"$\tau$ "
        self._xLabel += normalization + units
        self._yLabel = r"$\widetilde{n}$ " +\
                       yNormalization +\
                       yUnits

        # Set title
        if self._avg:
            title = r"$\mathrm{{Average \;{}}}$".format(blobOrHole[:-1])
        else:
            title = r"$\mathrm{{{}}}$".format(blobOrHole[:-1].capitalize())

        return title
    #}}}

    #{{{plotSaveShowTimeTrace
    def plotSaveShowTimeTrace(self):
        """
        Performs the plotting by calling the plot worker.
        """

        # Create the plot
        fig, ax = plt.subplots(figsize = SizeMaker.standard(w=4.0, a=0.5))

        self._plotWorker(ax, self._var, self._time, self._title)

        # Show and save
        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}

    #{{{_plotWorker
    def _plotWorker(self, ax, var, time, title, xbins=None):
        #{{{docstring
        """
        Performs axis plotting.

        Parameters
        ----------
        ax : axis
            The axis to plot on.
        var : array 1-d
            Variable to plot against time.
        time : array 1-d
            The time.
        title : str
            Title to  be used in the axis.
        xbins : [None|int]
            Number of bins to use on the x-axis
        """
        #}}}

        # Do the plotting
        ax.plot(time, var)
        ax.set_xlabel(self._xLabel)
        ax.set_ylabel(self._yLabel)

        # Set the title
        ax.set_title(title)

        # Make the plot look nice
        self._ph.makePlotPretty(ax, rotation = 45, legend = False, xbins=xbins)
    #}}}
#}}}

#{{{PlotBlobAndHoleTimeTraceDouble
class PlotBlobAndHoleTimeTraceDouble(PlotBlobOrHoleTimeTraceSingle):
    """
    Class which contains the time traces together with the plotting
    configuration of blobs and holes.
    """

    #{{{setData
    def setData(self, blobDict, holeDict):
        #{{{docstring
        """
        Sets the waiting time and pulse widths to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        blobDict : dict
            Dictionary of the blob.
            Contains the keys "n", "time" and "pos".
        holeDict : dict
            Dictionary of the hole.
            Contains the keys "n", "time" and "pos".
        """
        #}}}

        # Set the member data
        self._blobVar  = blobDict.pop("n")
        self._holeVar  = holeDict.pop("n")
        self._time     = blobDict.pop("time")
        self._pos      = blobDict.pop("pos")

        self._avg = True if self._time[0] < 0 else False

        avg = "-avg" if self._avg else ""

        fileName = "timeTrace-blobAndHole{}".format(avg)

        self._blobTitle = self._prepareLabels("blobs")
        self._holeTitle = self._prepareLabels("holes")

        self._setFileName(fileName)
    #}}}

    #{{{plotSaveShowTimeTrace
    def plotSaveShowTimeTrace(self):
        """
        Performs the plotting by calling the plot worker.
        """

        # Create the plot
        fig, (blobAx, holeAx) =\
                plt.subplots(ncols=2, figsize = SizeMaker.array(2,1))

        xbins = 5
        self._plotWorker(blobAx,\
                         self._blobVar,\
                         self._time,\
                         self._blobTitle,\
                         xbins = xbins,\
                         )
        self._plotWorker(holeAx,\
                         self._holeVar,\
                         self._time,\
                         self._holeTitle,\
                         xbins = xbins,\
                         )
        fig.subplots_adjust(wspace=0.75)

        # Show and save
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
    def __init__(self, *args, bins = "sqrt", **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets member data

        Parameters
        ----------
        bins : [int|str]
            The binning, see numpy.histogram for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._bins = bins
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
        figSize = SizeMaker.standard(w=6.0, a=1/3)
        fig, (wAx, pAx) =\
                plt.subplots(nrows=1, ncols=2, figsize = figSize)

        kwargs = {"normed":self._normed, "ec":"k", "fc":"g", "alpha":0.75}
        # Waiting time
        wAx.hist(self._pulseWidths, self._bins, **kwargs)
        wAx.set_xlabel(self._xLabel)
        wAx.set_ylabel(self._yLabelW)
        wAx.set_title(self._wTitle)

        # Pulse width
        pAx.hist(self._waitingTimes, self._bins, **kwargs)
        pAx.set_xlabel(self._xLabel)
        pAx.set_ylabel(self._yLabelP)
        pAx.set_title(self._pTitle)

        # Set the title
        fig.suptitle(self._title)

        # Modify title position
        t = fig.texts[0]
        pos = list(t.get_position())
        pos[1] = 1.0
        t.set_position(pos)
        t.set_va("bottom")

        # Make the plot look nice
        self._ph.makePlotPretty(wAx, rotation = 45, legend = False, xbins=7)
        self._ph.makePlotPretty(pAx, rotation = 45, legend = False, xbins=7)

        # Only integers on the yaxis
        wAx.get_yaxis().set_major_locator(MaxNLocator(integer=True))
        pAx.get_yaxis().set_major_locator(MaxNLocator(integer=True))

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.2, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
