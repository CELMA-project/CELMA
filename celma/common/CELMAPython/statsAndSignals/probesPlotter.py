#!/usr/bin/env python

""" Collection of plotting results for the probes """

from ..plotHelpers import plotNumberFormatter, seqCMap, seqCMap2
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.gridspec import GridSpec
import numpy as np
import os

# TODO: Add physical parameters calculation

#{{{PlotProbes
class PlotProbes(object):
    """
    Class which contains the probe data and the plotting configuration.

    NOTE: The probe position of the probes is only allowed to vary in rho
    """

    #{{{__init___
    def __init__(self             ,\
                 probes           ,\
                 showPlot  = False,\
                 savePlot  = False,\
                 extension = "png",\
                 savePath  = "."  ,\
                 pltSize   = None ,\
                 ):
        #{{{docstring
        """
        The constructor for the PlotProbes object.

        Sets the member data.

        NOTE: The probe position of the probes is only allowed to vary in rho

        Parameters
        ----------
        probes : Probes object
            Contains the data to plot.
        showPlot : bool
            If the plots should be displayed.
        savePlot : bool
            If the plots should be saved.
        extension : str
            Extension to use on the plots
        savePath : str
            Path to save destination. Must exist.
        pltSize : tuple
            Size of the plots given as (x, y)
        """
        #}}}

        # Set the member data
        self._probes    = probes
        self._showPlot  = showPlot
        self._savePlot  = savePlot
        self._extension = extension
        self._savePath  = savePath

        # Get the colors
        self._colors = seqCMap(np.linspace(0, 1, len(self._probes.probesKeys)))

        # Get the default title
        # NOTE: Theta should be the same for all probes
        #       (assume only rho changing)
        theta = self._probes.theta[self._probes.probesKeys[0]]*(180/np.pi)
        # Make from radians to degrees
        theta = plotNumberFormatter(theta, None)
        # NOTE: z should be the same for all probes (assume only rho changing)
        z     = plotNumberFormatter(self._probes.z[self._probes.probesKeys[0]], None)
        self._defaultTitle =\
            r"$\theta=${}$^{{\circ}},$ $z=${}".format(theta, z)

        # Set the plot size
        self._pltSize = pltSize
    #}}}

    #{{{plotTimeTrace
    def plotTimeTrace(self):
        """ Plots the time traces of the fluctuations."""

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        for nr, key in enumerate(self._probes.probesKeys):
            # Make the labels
            rho  = plotNumberFormatter(self._probes.rho[key], None)
            mean = plotNumberFormatter(self._probes.results[key]["mean"], None)
            var  = plotNumberFormatter(self._probes.results[key]["var"], None)
            label = r"$\rho=${:8} $\mu_n$={:8} $\sigma^2_n$={}".\
                    format(rho, mean, var)

            ax.plot(self._probes.time,\
                    self._probes.timeTraceOfVarFluct[key],\
                    color=self._colors[nr],\
                    label=label)

        # Set axis label
        ax.set_xlabel(r"$t\omega_{ci}$")
        ax.set_ylabel(r"$\tilde{n}$")

        fig.suptitle(self._defaultTitle)

        # Make the plot look nice
        self._makePlotPretty(ax)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._saveThePlot(fig, "timeTraces")

        plt.close(fig)
    #}}}

    #{{{plotPDFs
    def plotPDFs(self):
        """
        Plots the PDFs at all positions in the Probes object.
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        for nr, key in enumerate(self._probes.probesKeys):
            # Make the labels
            rho  = plotNumberFormatter(self._probes.rho[key], None)
            skew = plotNumberFormatter(self._probes.results[key]["skew"], None)
            kurt = plotNumberFormatter(self._probes.results[key]["kurtosis"], None)
            label = r"$\rho=${:8} $S_n$={:8} $K_n$={}".\
                    format(rho, skew, kurt)

            ax.plot(self._probes.results[key]["pdfX"],\
                    self._probes.results[key]["pdfY"],\
                    color=self._colors[nr],\
                    label=label)

        # Set logscale
        ax.set_yscale("log")

        # Set axis label
        ax.set_xlabel(r"$\tilde{n}$")
        ax.set_ylabel(r"$\mathrm{PDF}(\tilde{n})$")

        fig.suptitle(self._defaultTitle)

        # Make the plot look nice
        self._makePlotPretty(ax)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._saveThePlot(fig, "PDFs")

        plt.close(fig)
    #}}}

    #{{{plotAvgFluxThrougVolumeElement
    def plotAvgFluxThrougVolumeElement(self, uName, labelName):
        """
        Plots the average flux through a volume element

        uName : str
            Name of the velocity
        labelName : str
            Name of the velocity in a LaTeX plottable format (not
            including $)
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        gs      = GridSpec(nrows=3, ncols=1)
        totalAx = fig.add_subplot(gs[0])
        avgAx   = fig.add_subplot(gs[1], sharex=totalAx)
        fluctAx = fig.add_subplot(gs[2], sharex=totalAx)

        for nr, key in enumerate(self._probes.probesKeys):
            # Make the labels
            rho  = plotNumberFormatter(self._probes.rho[key], None)
            label = r"$\rho=${:8}".format(rho)

            totalAx.plot(self._probes.time,\
                         self._probes.results[key]["varAvgFlux" +\
                                                   uName.capitalize()],\
                         color=self._colors[nr],\
                         label=label)

            avgAx.  plot(self._probes.time,\
                         self._probes.results[key]["varAvgFluxAvg" +\
                                                   uName.capitalize()],\
                         color=self._colors[nr],\
                         label=label)

            fluctAx.plot(self._probes.time,\
                         self._probes.results[key]["varAvgFluxFluct" +\
                                                   uName.capitalize()],\
                         color=self._colors[nr],\
                         label=label)

        # Set axis label
        totalAx.set_ylabel(r"$n{}$".format(labelName))
        totalAx.tick_params(labelbottom="off")
        avgAx  .set_ylabel(r"$\bar{{n}}\bar{}$".format(labelName))
        avgAx  .tick_params(labelbottom="off")
        fluctAx.set_xlabel(r"$t\omega_{ci}$")
        fluctAx.set_ylabel(\
            r"$\langle\tilde{{n}}\tilde{}\rangle$".format(labelName))

        fig.suptitle(self._defaultTitle)

        # Make the plot look nice
        # Plot the legend
        leg = totalAx.legend(loc="best", fancybox = True, numpoints=1)
        leg.get_frame().set_alpha(0.5)

        axes = [totalAx, avgAx, fluctAx]

        for ax in axes:
            self._makePlotPretty(ax, prune = "both", rotation = 45)

        # Make sure no collision between the ticks
        fluctAx.xaxis.set_major_locator(MaxNLocator(prune="lower"))

        # Adjust the subplots
        fig.subplots_adjust(hspace=0)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._saveThePlot(fig, "flux{}".format(uName))

        plt.close(fig)
    #}}}

    #{{{plotZFFT
    def plotZFFT(self, positionKey, maxMode):
        """
        Plots the fourier transform at the given position.

        Parameters
        ----------
        positionKey : str
            The key to the self._probes.results which we would like to plot for.
            The key is on the form "dataXInd,dataYInd,dataZInd", where the
            indices refer to the keys in the simulated data.
        maxMode : int
            Number of modes to plot.
        """

        # Create figure and colors
        # NOTE: We use a different colormap here to distinguish from the
        #       positions
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)
        colors = seqCMap2(np.linspace(0, 1, maxMode))

        # Skip the offset mode
        for modeNr in range(1, maxMode+1):
            ax.plot(self._probes.time,\
                    self._probes.results[positionKey]["zFFT"][:, modeNr],\
                    color=colors[modeNr-1],\
                    label=r"$k_\theta={}$".format(modeNr))

        # Set logscale
        ax.set_yscale("log")

        # Set axis label
        ax.set_xlabel("Time []")
        ax.set_ylabel("Amplitude")

        rho = plotNumberFormatter(self._probes.rho[positionKey], None)
        z   = plotNumberFormatter(self._probes.z[positionKey], None)
        title  = r"$\rho=${}$,$ $z=${}".format(rho, z)

        fig.suptitle(title)

        # Make the plot look nice
        self._makePlotPretty(ax)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._saveThePlot(fig, "zFFT_at_{}".format(positionKey))

        plt.close(fig)
    #}}}

    #{{{_makePlotPretty
    def _makePlotPretty(self, ax, prune = "lower", rotation = "horizontal"):
        """
        Routine that fixes some beauty-mistakes in matplotlib

        Parameters
        ----------
        ax : axis
            The axis to fix.
        prune : str
            What ticks should be pruned.
        rotation : [str | int]
            Rotation of the x axis.
        """

        # Avoid silly top value (only for non-log axes)
        try:
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
        except:
            pass
        # Format the tick labels
        ax.get_xaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
        ax.get_yaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=rotation)
        # Plot the legend
        leg = ax.legend(loc="best", fancybox = True, numpoints=1)
        leg.get_frame().set_alpha(0.5)
        # Plot the grid
        ax.grid()
        # Make sure no collision between the ticks
        ax.xaxis.set_major_locator(MaxNLocator(prune=prune))
    #}}}

    #{{{_saveThePlot
    def _saveThePlot(self, fig, name):
        """
        Saves the figure

        Parameters
        ----------
        fig: figure
            The figure.
        name : str
            The name of the plot.
        """

        fileName = "{}.{}".\
            format(os.path.join(self._savePath, name), self._extension)

        fig.savefig(fileName,\
                    transparent = True    ,\
                    bbox_inches = "tight" ,\
                    pad_inches  = 0       ,\
                    )

        print("Saved to {}".format(fileName))
    #}}}
#}}}
