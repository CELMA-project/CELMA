#!/usr/bin/env python

""" Collection of plotting results for the probes """

from ..plotHelpers import plotNumberFormatter, seqCMap, seqCMap2
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
import numpy as np
import os

#{{{PlotProbes
class PlotProbes(object):
    """
    Class which contains the probe data and the plotting configuration.

    NOTE: The probe position of the probes is only allowed to vary in rho
    """

    #{{{__init___
    def __init__(self                     ,\
                 probes                   ,\
                 showPlot          = False,\
                 savePlot          = False,\
                 extension         = "png",\
                 savePath          = "."  ,\
                 pltSize           = None ,\
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
        self._colors = seqCMap(np.linspace(0, 1, len(probes.probesKeys)))

        # Make the default title
        # NOTE: Theta should be the same for all probes
        #       (assume only rho changing)
        theta = int(probes.theta[probes.probesKeys[0]])
        # NOTE: z should be the same for all probes (assume only rho changing)
        probes.helper.zTxtDict["value"] =\
                plotNumberFormatter(probes.z[probes.probesKeys[0]], None)
        self._defaultTitle = r"{}$,$ {}".\
             format(probes.helper.thetaTxtDict["constThetaTxt"].format(theta),\
                    probes.helper.zTxtDict["constZTxt"].format(probes.helper.zTxtDict))

        # Set default time label
        self._timeLabel = probes.helper.tTxtDict["tTxtLabel"].\
                          format(probes.helper.tTxtDict)

        # Set the variable label
        if probes.helper.convertToPhysical:
            self._varLabel = r"${}$ $[{}]$".\
                                  format(probes.varName, probes.varUnits)
            self._varLabelFLuct = r"$\tilde{{{}}}$ $[{}]$".\
                                  format(probes.varName, probes.varUnits)
        else:
            self._varLabel = r"${}{}$".\
                                  format(probes.varName, probes.varNormalization)
            self._varLabelFLuct = r"$\tilde{{{}}}{}$".\
                                  format(probes.varName, probes.varNormalization)

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
            self._probes.helper.rhoTxtDict["value"] =\
                    plotNumberFormatter(self._probes.rho[key], None)
            mean = plotNumberFormatter(self._probes.results[key]["mean"],\
                                       None,\
                                       precision=3)
            var  = plotNumberFormatter(self._probes.results[key]["var"],\
                                       None,\
                                       precision=3)
            rho = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                    format(self._probes.helper.rhoTxtDict)
            label = r"{}  $\mu_n$={:8} $\sigma^2_n$={}".\
                    format(rho, mean, var)

            ax.plot(self._probes.fluctTime,\
                    self._probes.timeTraceOfVarFluct[key],\
                    color=self._colors[nr],\
                    label=label)

        # Set axis labels
        ax.set_xlabel(self._timeLabel)
        ax.set_ylabel(self._varLabelFLuct)

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._probes.helper.makePlotPretty(ax, rotation = 45)

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
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "timeTraces"),\
                       self._extension)
            self._probes.helper.savePlot(fig, fileName, (self._leg,))

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
            self._probes.helper.rhoTxtDict["value"] =\
                    plotNumberFormatter(self._probes.rho[key], None)
            skew =\
                plotNumberFormatter(self._probes.results[key]["skew"], None)
            kurt =\
                plotNumberFormatter(self._probes.results[key]["kurtosis"], None)

            rho = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                    format(self._probes.helper.rhoTxtDict)
            label = r"{}   $S_n$={:8} $K_n$={}".\
                    format(rho, skew, kurt)

            ax.plot(self._probes.results[key]["pdfX"],\
                    self._probes.results[key]["pdfY"],\
                    color=self._colors[nr],\
                    label=label)

        # Set logscale
        ax.set_yscale("log")

        # Set axis label
        ax.set_xlabel(self._varLabelFLuct)
        ax.set_ylabel(r"$\mathrm{{PDF}}(\tilde{{{}}})$".\
                format(self._probes.varName))

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._probes.helper.makePlotPretty(ax)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "PDFs"),\
                       self._extension)
            self._probes.helper.savePlot(fig, fileName)

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
            self._probes.helper.rhoTxtDict["value"] =\
                    plotNumberFormatter(self._probes.rho[key], None)
            label = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                        format(self._probes.helper.rhoTxtDict)

            ax.plot(self._probes.results[key]["psdX"],\
                    self._probes.results[key]["psdY"],\
                    color=self._colors[nr],\
                    label=label)

        # Set logscale
        ax.set_yscale("log")

        # Set axis label
        if self._probes.helper.convertToPhysical:
            inverse = "$/Hz$"
            xlabel = "$\mathrm{f}$ $\mathrm{[Hz]}$"
        else:
            inverse = "$/${}".format(self._timeLabel)
            xlabel = "$(1/${}$)$".format(self._timeLabel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"${}^2${}".format(self._probes.varName, inverse))

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._probes.helper.makePlotPretty(ax, rotation = 45)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "PSDs"),\
                       self._extension)
            self._probes.helper.savePlot(fig, fileName)

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
            self._probes.helper.rhoTxtDict["value"] =\
                    plotNumberFormatter(self._probes.rho[key], None)
            label = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                        format(self._probes.helper.rhoTxtDict)

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
        if self._probes.helper.convertToPhysical:
            labelEnd = "[{}\mathrm{{ms}}^{{-1}}]".format(self._probes.varUnits)
        else:
            labelEnd = "{}c_s".format(self._probes.varNormalization)

        if pltFluct:
            axes["fluctAx"].set_ylabel(\
                r"$\langle\tilde{{{}}}\tilde{{{}}}\rangle{}$".\
                    format(self._probes.varName, labelName, labelEnd))
        if pltTotal:
            axes["totalAx"].set_ylabel(r"${}{}{}$".\
                    format(self._probes.varName, labelName, labelEnd))
        if pltAvg:
            axes["avgAx"].set_ylabel(r"$\bar{{{}}}\bar{}{}$".\
                    format(self._probes.varName, labelName, labelEnd))

        pltAxes[0].set_title(self._defaultTitle)

        # Make the plot look nice
        # Plot the legend
        leg = pltAxes[0].legend(loc="best", fancybox = True, numpoints=1)
        leg.get_frame().set_alpha(0.5)

        for ax in pltAxes:
            self._probes.helper.makePlotPretty(ax, yprune = "both", rotation = 45)

        for ax in pltAxes[0:-1]:
            avgAx.tick_params(labelbottom="off")

        # Make sure no collision between the ticks
        pltAxes[-1].xaxis.set_major_locator(MaxNLocator(prune="lower"))
        pltAxes[-1].set_xlabel(self._timeLabel)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0)

        # Manual tweeking as we want legends outside the plot
        fig.tight_layout(rect=[0,0,0.9,0.95])


        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "flux{}".format(uName)),\
                       self._extension)
            self._probes.helper.savePlot(fig, fileName)

        plt.close(fig)
    #}}}

    #{{{plotZFFT
    def plotZFFT(self, positionKey, maxMode):
        r"""
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
        ax.set_xlabel(self._timeLabel)
        ax.set_ylabel(self._varLabel)

        self._probes.helper.rhoTxtDict["value"] =\
                plotNumberFormatter(self._probes.rho[positionKey], None)
        rho = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                    format(self._probes.helper.rhoTxtDict)
        self._probes.helper.zTxtDict  ["value"] =\
                plotNumberFormatter(self._probes.z[positionKey], None)
        z   = self._probes.helper.zTxtDict["constZTxt"].\
                    format(self._probes.helper.zTxtDict)
        title  = r"{}$,$ {}".format(rho, z)

        ax.set_title(title)

        # Make the plot look nice
        self._probes.helper.makePlotPretty(ax, rotation = 45)
        fig.tight_layout()

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath,\
                       "zFFT_at_{}".format(positionKey)),\
                       self._extension)
            self._probes.helper.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
