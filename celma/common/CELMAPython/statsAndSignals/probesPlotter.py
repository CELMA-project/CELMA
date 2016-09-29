#!/usr/bin/env python

""" Collection of plotting results for the probes """

from ..plotHelpers import plotNumberFormatter, seqCMap, seqCMap2
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
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
        self._alpha     = 0.7

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
            fluctMean =\
                plotNumberFormatter(self._probes.results[key]["fluctMean"],\
                                    None,\
                                    precision=3)
            fluctVar =\
                plotNumberFormatter(self._probes.results[key]["fluctVar"],\
                                    None,\
                                    precision=3)
            rho = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                    format(self._probes.helper.rhoTxtDict)
            label = (r"{}  $ \mu_{{\tilde{{n}}}}$={:8} "
                     r"$\sigma^2_{{\tilde{{n}}}}$={}").\
                             format(rho, fluctMean, fluctVar)

            ax.plot(self._probes.fluctTime,\
                    self._probes.\
                        timeTraceOfVar[key][self._probes.tIndSaturatedTurb:],\
                    color=self._colors[nr],\
                    label=label)

        # Set axis labels
        ax.set_xlabel(self._timeLabel)
        ax.set_ylabel(self._varLabel)

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
            fluctSkew =\
               plotNumberFormatter(self._probes.results[key]["fluctSkew"], None)
            fluctKurt =\
               plotNumberFormatter(self._probes.results[key]["fluctKurt"], None)

            rho = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                    format(self._probes.helper.rhoTxtDict)
            label = r"{}   $S_{{\tilde{{n}}}}$={:8} $K_{{\tilde{{n}}}}$={}".\
                    format(rho, fluctSkew, fluctKurt)

            ax.plot(self._probes.results[key]["pdfX"],\
                    self._probes.results[key]["pdfY"],\
                    color=self._colors[nr],\
                    label=label,\
                    alpha=self._alpha)

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
            # Clip the very first point as this is rediculously low
            ax.plot(self._probes.results[key]["psdX"][1:],\
                    self._probes.results[key]["psdY"][1:],\
                    color=self._colors[nr],\
                    label=label,\
                    alpha=self._alpha)

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

    #{{{plotAvgFluxThroughVolumeElement
    def plotAvgFluxThroughVolumeElement(self,\
                                        uName,\
                                        labelName,\
                                       ):
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

        totalPlots = 0

        gs = GridSpec(nrows=len(self._probes.probesKeys), ncols=1)

        axes = []
        for nr, key in enumerate(self._probes.probesKeys):
            if nr == 0:
                axes.append(fig.add_subplot(gs[nr]))
            else:
                axes.append(fig.add_subplot(gs[nr], sharex=axes[0]))

            # Make the labels
            self._probes.helper.rhoTxtDict["value"] =\
                    plotNumberFormatter(self._probes.rho[key], None)
            label = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                        format(self._probes.helper.rhoTxtDict)

            axes[nr].plot(self._probes.fluctTime,\
                          self._probes.\
                             results[key]["varAvgFluxFluct" +\
                                         uName.capitalize()],\
                          color=self._colors[nr],\
                          label=label)
            # Plot the legends
            leg = axes[nr].legend(loc="upper right", fancybox = True, numpoints=1)
            leg.get_frame().set_alpha(0.5)

        # Set axis label
        if self._probes.helper.convertToPhysical:
            labelEnd =\
                "[{}\mathrm{{ms}}^{{-1}}]".format(self._probes.varUnits)
        else:
            labelEnd =\
                "{}c_s".format(self._probes.varNormalization)

        midAx = round(len(axes)/2)
        axes[midAx].set_ylabel(\
            r"$\langle\tilde{{{}}}\tilde{{{}}}\rangle{}$".\
                format(self._probes.varName, labelName, labelEnd))

        axes[0].set_title(self._defaultTitle)

        # Make the plot look nice
        for ax in axes:
            self._probes.helper.makePlotPretty(ax, yprune = "both",\
                                               rotation = 45)

        for ax in axes[0:-1]:
            ax.tick_params(labelbottom="off")

        # Make sure no collision between the ticks
        axes[-1].xaxis.set_major_locator(MaxNLocator(prune="lower"))
        axes[-1].set_xlabel(self._timeLabel)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0)

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
    def plotZFFT(self, positionKey, maxMode, clip=True, plotLinear=True):
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
        clip : bool
            If True, the first three points (typically where the
            perturbed noise decays) is removed from the plot.
        plotLinear : bool
            If True, the linear phase will be plotted as well
        """

        #{{{plotFunc
        def plotFunc(linear=False):
            """
            The function which plots and saves the FFT plot

            Parameters
            ----------
            linear : bool
                If True, "linear" will be appended to the plot name.
            """
            # Create figure and colors
            # NOTE: We use a different colormap here to distinguish from the
            #       positions
            fig = plt.figure(figsize = self._pltSize)
            ax  = fig.add_subplot(111)
            colors = seqCMap2(np.linspace(0, 1, maxMode))

            # Skip the offset mode
            for modeNr in range(1, maxMode+1):
                # Clip where the modes has been added
                clip = 3
                if linear:
                    linearClip = self._probes.results[positionKey]["zFFTLinearIndex"]
                else:
                    linearClip = None

                ax.plot(self._probes.time[clip:linearClip],\
                        self._probes.results[positionKey]["zFFT"]\
                            [clip:linearClip, modeNr],\
                        color=colors[modeNr-1],\
                        label=r"$k_\theta={}$".format(modeNr),
                        alpha=self._alpha)

            # Set logscale
            try:
                ax.set_yscale("log")
            except ValueError as er:
                if "no positive values" in er.args[0]:
                    message = ("{0}{1}WARNING: Plot of FFT failed as no "
                               "positive values{1}{0}".format("\n"*2, "!"*3))
                    print(message)
                    return


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
                if linear:
                    linearOrNot = "Linear"
                else:
                    linearOrNot = ""
                fileName = "{}{}.{}".\
                    format(os.path.join(self._savePath,\
                                "zFFT_at_{}".format(positionKey)),\
                           linearOrNot,\
                           self._extension)
                self._probes.helper.savePlot(fig, fileName)

            plt.close(fig)
        #}}}

        # Call the plot function
        plotFunc(linear=False)

        if plotLinear:
            plotFunc(linear=True)
    #}}}
#}}}
