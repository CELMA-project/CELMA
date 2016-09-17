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
                 convertToPhysical   =False,\
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
        convertToPhysical : str
            Whether normalized or physical units should be used.
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

        # Get labels
        # String formatting
        self._tTxtDict   = {"normalization":tNormalization  , "units":tUnits}
        self._rhoTxtDict = {"normalization":rhoNormalization, "units":rhoUnits}
        self._zTxtDict   = {"normalization":zNormalization  , "units":zUnits}

        self._rhoTxt   = r"$\rho{0[normalization]}$".format(self._rhoTxtDict)
        self._thetaTxt = r"$\theta={:d}^{{\circ}}$"
        self._zTxt     = r"$z{0[normalization]}$".format(self._zTxtDict)

        # Expand the dictionaries
        self._rhoTxtDict['rhoTxt'] = self._rhoTxt
        self._zTxtDict['zTxt']     = self._zTxt

        if self.convertToPhysical:
            self._rhoTxtLabel = "{0[rhoTxt]} $[{0[units]}]$".\
                    format(self._rhoTxtDict)
            self._zTxtLabel   = "{0[zTxt]} $[{0[units]}]$".\
                    format(self._zTxtDict)

            self._constRhoTxt = r"{0[rhoTxt]} $=$ {0[value]} ${0[units]}$"
            self._constZTxt   = r"{0[zTxt]} $=$ {0[value]} ${0[units]}$"
            self._tTxt        =\
                r"$\mathrm{{t}}{0[normalization]}$ $=$ {0[value]} ${0[units]}$"
        else:
            self._rhoTxtLabel = "{0[rhoTxt]}".format(self._rhoTxtDict)
            self._zTxtLabel   = "{0[zTxt]}"  .format(self._zTxtDict)

            self._constRhoTxt = r"{0[rhoTxt]} $=$ {0[value]}"
            self._constZTxt   = r"{0[zTxt]} $=$ {0[value]}"
            self._tTxt        =\
                r"$t{0[normalization]}$ $=$ {0[value]}"

        # Set the plot size
        self._pltSize = pltSize
    #}}}

    #{{{plotTimeTrace
    def plotTimeTrace(self):
        """ Plots the time traces of the fluctuations."""

        # Manual tweeking as we want legends outside the plot
        pltSize = (self._pltSize[0] + 10, self._pltSize[1])
        # Create the plot
        fig = plt.figure(figsize = pltSize)
        ax  = fig.add_subplot(111)

        for nr, key in enumerate(self._probes.probesKeys):
            # Make the labels
            rho  = plotNumberFormatter(self._probes.rho[key], None)
            mean = plotNumberFormatter(self._probes.results[key]["mean"],\
                                       None,\
                                       precision=3)
            var  = plotNumberFormatter(self._probes.results[key]["var"],\
                                       None,\
                                       precision=3)
            label = r"$\rho=${:8} $\mu_n$={:8} $\sigma^2_n$={}".\
                    format(rho, mean, var)

            ax.plot(self._probes.fluctTime,\
                    self._probes.timeTraceOfVarFluct[key],\
                    color=self._colors[nr],\
                    label=label)

        # Set axis label
        # FIXME: Units in the old fashioned way
        ax.set_xlabel(r"$t\omega_{ci}$")
        ax.set_ylabel(r"$\tilde{n}$")

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._makePlotPretty(ax, rotation = 45)

        # Manual tweeking as we want legends outside the plot
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        self._leg = ax.legend(loc="center left",\
                              fancybox = True,\
                              numpoints=1,\
                              bbox_to_anchor=(1, 0.5),\
                              )
        self._leg.get_frame().set_alpha(0.5)
        # Manual tweeking as we want legends outside the plot
        fig.tight_layout(rect=[0,0,0.7,1])

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

        NOTE: The PDF is a probability, and thus has no units.
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
        # FIXME: Units of the variables here
        ax.set_xlabel(r"$\tilde{n}$")
        # FIXME: Units in the old fashioned way
        ax.set_ylabel(r"$\mathrm{PDF}(\tilde{n})$")

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._makePlotPretty(ax)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._saveThePlot(fig, "PDFs")

        plt.close(fig)
    #}}}

    #{{{plotPSDs
    def plotPSDs(self):
        """
        Plots the PSDs at all positions in the Probes object.

        NOTE: The units will be in variableUnits**2/Hz for
              non-normalized variables.
              The normalization would be
              variableNormalization**2/(tOmegaCI)
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        for nr, key in enumerate(self._probes.probesKeys):
            # Make the labels
            rho  = plotNumberFormatter(self._probes.rho[key], None)
            label = r"$\rho=${}".format(rho)

            ax.plot(self._probes.results[key]["psdX"],\
                    self._probes.results[key]["psdY"],\
                    color=self._colors[nr],\
                    label=label)

        # Set logscale
        ax.set_yscale("log")

        # Set axis label
        # FIXME: Herz or so
        ax.set_xlabel(r"${}$")
        # FIXME: Units as described above
        ax.set_ylabel(r"${}^2/{}$")

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._makePlotPretty(ax)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._saveThePlot(fig, "PSDs")

        plt.close(fig)
    #}}}

    #{{{plotAvgFluxThrougVolumeElement
    def plotAvgFluxThrougVolumeElement(self,\
                                       uName,\
                                       labelName,\
                                       pltFluct = True,\
                                       pltTotal = False,\
                                       pltAvg   = False,\
                                       ):
        """
        Plots the average flux through a volume element

        uName : str
            Name of the velocity
        labelName : str
            Name of the velocity in a LaTeX plottable format (not
            including $)
        pltFluct : bool
            If the average flux fluctuation is to be plotted.
        pltTotal : bool
            If the total average flux is to be plotted.
        pltAvg : bool
            If the average flux average is to be plotted.
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)

        totalPlots = 0
        if pltFluct:
            totalPlots += 1
        if pltTotal:
            totalPlots += 1
        if pltAvg:
            totalPlots += 1

        gs = GridSpec(nrows=totalPlots, ncols=1)

        axes = {}
        pltAxes = []
        curPlot = 0
        if pltFluct:
            axes["fluctAx"] = fig.add_subplot(gs[curPlot])
            curPlot += 1
            pltAxes.append(axes["fluctAx"])
        if pltTotal:
            if pltFluct:
                axes["totalAx"] =\
                        fig.add_subplot(gs[curPlot], sharex=axes["fluctAx"])
            else:
                axes["totalAx"] = fig.add_subplot(gs[curPlot])
            curPlot += 1
            pltTotal.append(axes["totalAx"])
        if pltAvg:
            if pltFluct:
                axes["avgAx"] =\
                        fig.add_subplot(gs[curPlot], sharex=axes["fluctAx"])
            elif pltTotal:
                axes["avgAx"] =\
                        fig.add_subplot(gs[curPlot], sharex=axes["totalAx"])
            else:
                axes["avgAx"] = fig.add_subplot(gs[curPlot])

            pltTotal.append(axes["avgAx"])
        for nr, key in enumerate(self._probes.probesKeys):
            # Make the labels
            rho  = plotNumberFormatter(self._probes.rho[key], None)
            label = r"$\rho=${:8}".format(rho)

            if pltFluct:
                axes["fluctAx"].plot(self._probes.fluctTime,\
                                     self._probes.\
                                        results[key]["varAvgFluxFluct" +\
                                                    uName.capitalize()],\
                                     color=self._colors[nr],\
                                     label=label)
            if pltTotal:
                axes["totalAx"].plot(self._probes.time,\
                                     self._probes.\
                                        results[key]["varAvgFlux" +\
                                                    uName.capitalize()],\
                                     color=self._colors[nr],\
                                     label=label)
            if pltAvg:
                axes["avgAx"].plot(self._probes.fluctTime,\
                                   self._probes.\
                                    results[key]["varAvgFluxAvg" +\
                                                 uName.capitalize()],\
                                   color=self._colors[nr],\
                                   label=label)


        # Set axis label
        # FIXME: Units (need to multiply with m\s or cs)
        if pltFluct:
            axes["fluctAx"].set_ylabel(\
                r"$\langle\tilde{{n}}\tilde{}\rangle$".format(labelName))
        if pltTotal:
            axes["totalAx"].set_ylabel(r"$n{}$".format(labelName))
        if pltAvg:
            axes["avgAx"].set_ylabel(r"$\bar{{n}}\bar{}$".format(labelName))

        pltAxes[0].set_title(self._defaultTitle)

        # Make the plot look nice
        # Plot the legend
        leg = pltAxes[0].legend(loc="best", fancybox = True, numpoints=1)
        leg.get_frame().set_alpha(0.5)

        for ax in pltAxes:
            self._makePlotPretty(ax, prune = "both", rotation = 45)

        for ax in pltAxes[0:-1]:
            avgAx.tick_params(labelbottom="off")

        # Make sure no collision between the ticks
        pltAxes[-1].xaxis.set_major_locator(MaxNLocator(prune="lower"))
        # FIXME: Units and stuff
        pltAxes[-1].set_xlabel(r"$t\omega_{ci}$")

        # Adjust the subplots
        fig.subplots_adjust(hspace=0)

        # Manual tweeking as we want legends outside the plot
        fig.tight_layout(rect=[0,0,0.9,0.95])


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

        NOTE: As we are taking the fourier transform in the
              dimensionless theta direction, that means that the units
              will remain the same as

              \hat{f}(x) = \int_\infty^\infty f(x) \exp(-i2\pi x\xi) d\xi

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
        # FIXME: Units and stuff here
        ax.set_xlabel("$t\omega_{ci}$")
        # FIXME: Units and stuff here
        ax.set_ylabel("$\mathrm{Amplitude}$")

        rho = plotNumberFormatter(self._probes.rho[positionKey], None)
        z   = plotNumberFormatter(self._probes.z[positionKey], None)
        title  = r"$\rho=${}$,$ $z=${}".format(rho, z)

        ax.set_title(title)

        # Make the plot look nice
        self._makePlotPretty(ax, rotation = 45)
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
        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=rotation)
        ax.get_xaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
        ax.get_yaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
        # Plot the legend
        self._leg = ax.legend(loc="best",\
                              fancybox = True,\
                              numpoints=1,\
                              )
        self._leg.get_frame().set_alpha(0.5)
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
                    transparent = True             ,\
                    bbox_inches = "tight"          ,\
                    bbox_extra_artists=(self._leg,),\
                    pad_inches  = 0                ,\
                    )

        print("Saved to {}".format(fileName))
    #}}}
#}}}
