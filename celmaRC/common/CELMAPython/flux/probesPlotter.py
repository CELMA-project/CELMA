#!/usr/bin/env python

""" Collection of plotting results for the probes """

from ..plotHelpers import plotNumberFormatter, seqCMap2, seqCMap3
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
        self._alpha     = 0.7

        # Get the colors
        self._colors = seqCMap3(np.linspace(0, 1, len(probes.probesKeys)))

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
        fig.tight_layout(rect=(0,0,0.7,1))

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
            if np.isnan(self._probes.results[key]["pdfY"][0]):
                message = ("{0}{1} WARNING No PDF created for {2}{1}{0}")
                print(message.format("\n", "!"*3, key))
                continue
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

        if len(tuple(ax.get_lines())) == 0:
            message = "{0}{1}WARNING No PDFs to plot. Returning{1}{0}"
            print(message.format("\n", "!"*3))
            return

        # Set logscale
        ax.set_yscale("log")

        # Set axis label
        ax.set_xlabel(self._varLabelFLuct)
        if self._probes.helper.convertToPhysical:
            ax.set_ylabel(r"$\mathrm{{PDF}}(\tilde{{{}}}{})$".\
                format(self._probes.varName, self._probes.varNormalization))
        else:
            ax.set_ylabel(r"$\mathrm{{PDF}}(\tilde{{{}}})$".\
                    format(self._probes.varName))

        ax.set_title(self._defaultTitle)

        # Make the plot look nice
        self._probes.helper.makePlotPretty(ax, rotation = 45)
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

        #{{{plotFunc
        def plotFunc(mode = "normal"):
            if self._probes.results[self._probes.probesKeys[0]]["psdX"] is None:
                return

            # Create the plot
            fig = plt.figure(figsize = self._pltSize)
            ax  = fig.add_subplot(111)

            if mode == "normal":
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

            elif mode == "dB":
                # Find max:
                curMax = 0
                for nr, key in enumerate(self._probes.probesKeys):
                    if np.max(self._probes.results[key]["psdY"][1:]) > curMax:
                        curMax = np.max(self._probes.results[key]["psdY"][1:])

                # Make the plots
                for nr, key in enumerate(self._probes.probesKeys):
                    # Make the labels
                    self._probes.helper.rhoTxtDict["value"] =\
                            plotNumberFormatter(self._probes.rho[key], None)
                    label = self._probes.helper.rhoTxtDict["constRhoTxt"].\
                                format(self._probes.helper.rhoTxtDict)
                    # Clip the very first point as this is rediculously low
                    ax.plot(self._probes.results[key]["psdX"][1:],\
                            np.log10(\
                                self._probes.results[key]["psdY"][1:]/\
                                curMax),\
                            color=self._colors[nr],\
                            label=label,\
                            alpha=self._alpha)


            if mode == "normal":
                # Set logscale
                ax.set_yscale("log")
            else:
                ax.set_xscale("log")

            # Set axis label
            if self._probes.helper.convertToPhysical:
                inverse = "$/\mathrm{Hz}]$"
                xlabel = "$\mathrm{f}$ $\mathrm{[Hz]}$"
                if mode == "normal":
                    ax.set_ylabel(r"$\mathrm{{PSD}}$ $[({})^2${}".\
                            format(self._probes.varUnits, inverse))
                elif mode == "dB":
                    ax.set_ylabel(r"$\mathrm{dB}$")
            else:
                inverse = "{}".format(self._timeLabel)
                xlabel = "$(1/${}$)$".format(self._timeLabel)
                if mode == "normal":
                    ax.set_ylabel(r"${}{}^2${}".\
                            format(self._probes.varName,\
                                   self._probes.varNormalization,\
                                   inverse))
                elif mode == "dB":
                    ax.set_ylabel(r"$\mathrm{PSD}$ $\mathrm{dB}$")

            ax.set_xlabel(xlabel)

            ax.set_title(self._defaultTitle)

            # Make the plot look nice
            self._probes.helper.makePlotPretty(ax, rotation = 45)
            fig.tight_layout()

            if self._showPlot:
                plt.show()

            if self._savePlot:
                name = os.path.join(self._savePath, "PSDs")
                if mode == "dB":
                    name += "dB"
                fileName = "{}.{}".format(name, self._extension)
                self._probes.helper.savePlot(fig, fileName)

            plt.close(fig)
        #}}}

        plotFunc(mode = "normal")
        plotFunc(mode = "dB")
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
                                              rotation = 45, loc="upper right")

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
    def plotZFFT(self, positionKey, maxMode, clip=True,\
                 plotNonSaturated=True, plotLinear=True):
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
        plotNonSaturated : bool
            If True, the non saturated phase will be plotted as well
        plotLinear : bool
            If True, the linear phase will be plotted as well
        """

        #{{{plotFunc
        def plotFunc(plotRange="normal"):
            """
            The function which plots and saves the FFT plot

            Parameters
            ----------
            plotRange : ["normal"|"nonSaturated"|"linear"]
                * "normal" plots without restricting the end time
                * "nonSaturated" plots all the way until the saturated state
                * "linear" plots only the linear stage
            """
            # Create figure and colors
            # NOTE: We use a different colormap here to distinguish from the
            #       positions
            fig = plt.figure(figsize = self._pltSize)
            ax  = fig.add_subplot(111)
            colors = seqCMap2(np.linspace(0, 1, maxMode))

            # Number of points
            N = self._probes.results[positionKey]["zFFT"].shape[1]

            # Skip the offset mode
            for modeNr in range(1, maxMode+1):
                # Clip where the modes has been added
                clip = 3
                if plotRange == "normal":
                    endClip = None
                elif plotRange == "nonSaturated":
                    endClip = self._probes.results[positionKey]\
                                                  ["zFFTNonSaturatedIndex"]
                elif plotRange == "linear":
                    endClip = self._probes.results[positionKey]\
                                                  ["zFFTLinearIndex"]

                #{{{ NOTE: We are dealing with a real signal:
                #          As the fourier transform breaks the signal up
                #          in cisoids there will be one part of the
                #          signal in the positive rotating ciscoid and
                #          one in the negative (negative frequencies)
                #          for a given mode number. We need to take into
                #          account both in order to calculate the
                #          amplitude. As the signal is real only one of
                #          the phase sifts are needed. Notice that for a
                #          real signal the imaginary part occurs as a
                #          complex conjugate pair
                # http://dsp.stackexchange.com/questions/431/what-is-the-physical-significance-of-negative-frequencies?noredirect=1&lq=1
                # http://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
                #}}}
                # Magnitude of the signal
                # https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition
                posFreq = np.abs(self._probes.results[positionKey]["zFFT"]\
                            [clip:endClip, modeNr ])
                negFreq = np.abs(self._probes.results[positionKey]["zFFT"]\
                            [clip:endClip, -modeNr])
                modeMag = (posFreq + negFreq)/N

                ax.plot(self._probes.time[clip:endClip],\
                        modeMag,\
                        color=colors[modeNr-1],\
                        label=r"$m_\theta={}$".format(modeNr),
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

                if plotRange == "normal":
                    extraString = ""
                elif plotRange == "nonSaturated":
                    extraString = "NonSaturated"
                elif plotRange == "linear":
                    extraString = "Linear"

                fileName = "{}{}.{}".\
                    format(os.path.join(self._savePath,\
                                "zFFT_at_{}".format(positionKey)),\
                           extraString,\
                           self._extension)
                self._probes.helper.savePlot(fig, fileName)

            plt.close(fig)
        #}}}

        # Call the plot function
        plotFunc(plotRange="normal")
        if plotNonSaturated:
            plotFunc(plotRange="nonSaturated")
        if plotLinear:
            plotFunc(plotRange="linear")
    #}}}
#}}}
