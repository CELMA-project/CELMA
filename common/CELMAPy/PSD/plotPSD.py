#!/usr/bin/env python

"""Class for PSD plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import SizeMaker, plotNumberFormatter, seqCMap, seqCMap3
from matplotlib.ticker import FuncFormatter
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotPSD
class PlotPSD(PlotSuperClass):
    """
    Class which contains the PSD data and the plotting configuration.
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
    def setData(self, PSD, mode):
        #{{{docstring
        """
        Sets the PSD to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        PSD : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"pdfX":pdfX, "pdfY":"pdfY"}
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        """
        #}}}

        # Set the member data
        self._PSD = PSD
        self._mode = mode

        # Obtain the varname
        ind  = list(PSD.keys())[0]
        keys = PSD[ind].keys()
        self._varName = [var[:-4] for var in keys if "PSD" in var][0]

        # Obtain the color (pad away brigthest colors)
        pad = 3
        self._colors = seqCMap3(np.linspace(0, 1, len(PSD.keys())+pad))

        self._prepareLabels()

        # Set the labels
        pltVarName = self._ph.getVarPltName(self._varName)
        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, "PSD", self._fluctName))

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """
        # Set var label templates
        # NOTE: The units will be in variableUnits**2/Hz for
        #       non-normalized variables.
        #       The normalization would be
        #       variableNormalization**2/(tOmegaCI)
        #       However, we will use a dB approach, so all units cancel
        psdStr = r"\mathrm{{PSD}}("

        # Get normalOrFluct
        if self._mode == "normal":
            normalOrFluct = r"{{}}"
            self._fluctName = ""
        elif self._mode == "fluct":
            normalOrFluct = r"\widetilde{{{}}}"
            self._fluctName = "fluct"
        else:
            message = "'{}'-mode not implemented.".format(self._mode)
            raise NotImplementedError(message)

        #  Normalized and non-normalized labels are treated differently
        #  due to the paranthesis and the normalization
        units = "[\mathrm{{dB}}]"
        if self.uc.convertToPhysical:
            self._varLabelTemplate = r"$"+\
                                     psdStr +\
                                     normalOrFluct +\
                                     r")$" +\
                                     r" $" + units + r"$"
            self._freqLabel        = r"$[Hz]$"
        else:
            self._varLabelTemplate = r"$"+\
                                     psdStr +\
                                     normalOrFluct +\
                                     "{normalization}" +\
                                     r")$" +\
                                     r" $" + units + r"$"
            self._freqLabel        = r"$1/t \omega_{ci}$"
    #}}}

    #{{{setData2D
    def setData2D(self, PSD, mode, plotLimits):
        #{{{docstring
        """
        Sets the PSD to be plotted.

        This function also sets:
            * The variable labels
            * The colors
            * The contourplot keyword arguments
            * The save name

        Parameters
        ----------
        PSD : dict
            Dictionary with the keys
                * freqPosMatrix - array-2d on with dimension (rho, freq)
                * thetaPos      - integer with the fixed theta position
                * zPos          - integer with the fixed z position
                * FREQ          - meshgrid with the frequencies
                * RHO           - meshgrid with the rho coordinate
                * varName       - the variabel name
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        plotLimits : dict
            Dictionary on the form
            {"xlim":(min,max), "ylim":(min,max), "zlim":(min,max)}
            The dictionary values may be None rather than a tuple.
        """
        #}}}

        # Magic number
        nCont = 100

        # Set the member data
        self._PSD        = PSD
        self._mode       = mode
        self._plotLimits = plotLimits

        # Obtain the varname
        self._varName = self._PSD.pop("varName")

        # Set the z limits
        if self._plotLimits["zlim"] is None:
            vMin = np.min(self._PSD["freqPosMatrix"])
            vMax = np.max(self._PSD["freqPosMatrix"])
        else:
            vMin = self._plotLimits["zlim"][0]
            vMax = self._plotLimits["zlim"][1]

        levels = np.linspace(vMin, vMax, nCont, endpoint = True)

        # Make the contourf keyword arguments
        self._cfKwargs = {\
                          "vmax"   : vMax   ,\
                          "vmin"   : vMin   ,\
                          "levels" : levels ,\
                          "cmap"   : seqCMap,\
                          "extend" : "min"  ,\
                          "zorder" : -20    ,\
                         }

        self._prepareLabels()

        # Set the labels
        pltVarName = self._ph.getVarPltName(self._varName)
        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])

        # Set the x-axis label
        self._rhoLabel = self._ph.rhoTxtDict["rhoTxtLabel"]

        # Set the title
        # Set values
        z     = self._PSD.pop("zPos")
        theta = self._PSD.pop("thetaPos")
        self._ph.zTxtDict    ["value"] =\
                plotNumberFormatter(float(z), None)
        self._ph.thetaTxtDict["value"] =\
                plotNumberFormatter(float(theta), None)

        # Make the const values
        self._title = (r"{}$,$ {}").\
                        format(\
                            self._ph.thetaTxtDict["constThetaTxt"].\
                                format(self._ph.thetaTxtDict),\
                            self._ph.zTxtDict["constZTxt"].\
                                format(self._ph.zTxtDict),\
                              )

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, "PSD2D", self._fluctName))

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{plotSaveShowPSD
    def plotSaveShowPSD(self):
        """
        Plots the probability density function.

        NOTE: The desibel use here is dividing by its own max, but not
              the global max.
        """

        # Create the plot
        fig = plt.figure(figsize = SizeMaker.standard())
        ax  = fig.add_subplot(111)

        keys = sorted(self._PSD.keys())

        for key, color in zip(keys, self._colors):
            # Make the label
            rho, theta, z = key.split(",")

            # Set values
            self._ph.rhoTxtDict  ["value"] =\
                    plotNumberFormatter(float(rho), None)
            self._ph.zTxtDict    ["value"] =\
                    plotNumberFormatter(float(z), None)
            self._ph.thetaTxtDict["value"] =\
                    plotNumberFormatter(float(theta), None)

            # Make the const values
            label = (r"{}$,$ {}$,$ {}").\
                    format(\
                        self._ph.rhoTxtDict  ["constRhoTxt"].\
                            format(self._ph.rhoTxtDict),\
                        self._ph.thetaTxtDict["constThetaTxt"].\
                            format(self._ph.thetaTxtDict),\
                        self._ph.zTxtDict["constZTxt"].\
                            format(self._ph.zTxtDict),\
                          )

            # Clip the very first point as this is very low
            xVals = self._PSD[key]["{}PSDX".format(self._varName)][1:]
            yVals = np.log10(       self._PSD[key]["{}PSDY".\
                                     format(self._varName)][1:]/\
                             np.max(self._PSD[key]["{}PSDY".\
                                    format(self._varName)][1:]))

            # Plot
            ax.plot(xVals, yVals, color=color, label=label)

        # Use logarithmic scale
        ax.set_xscale("log")

        # Set axis labels
        ax.set_xlabel(self._freqLabel)
        ax.set_ylabel(self._varLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(ax, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "PSD"),\
                       self._extension)
            self._ph.savePlot(fig, fileName)

        plt.close(fig)
    #}}}

    #{{{plotSaveShowPSD2D
    def plotSaveShowPSD2D(self):
        """
        Plots the probability density function as a countourplot.
        """

        # Create the plot
        fig = plt.figure(figsize = SizeMaker.standard(w=6, a=0.3))
        ax  = fig.add_subplot(111)

        # Plot
        CP = ax.contourf(\
            self._PSD["RHO"], self._PSD["FREQ"], self._PSD["freqPosMatrix"],\
            **self._cfKwargs)

        if self._plotLimits["xlim"] is not None:
            ax.set_xlim(self._plotLimits["xlim"])
        if self._plotLimits["ylim"] is not None:
            ax.set_ylim(self._plotLimits["ylim"])

        cbar = fig.colorbar(CP, format = FuncFormatter(plotNumberFormatter))
        cbar.set_label(self._varLabel)
        # Out of range
        CP.set_clim(*self._plotLimits["zlim"])
        CP.cmap.set_under("k")

        # Set rasterization order
        ax.set_rasterization_zorder(-10)

        # Set axis labels
        ax.set_xlabel(self._rhoLabel)
        ax.set_ylabel(self._freqLabel)

        ax.set_title(self._title)

        # Make the plot look nice
        self._ph.makePlotPretty(ax,\
                                xprune   = "both",\
                                yprune   = "both",\
                                xbins    = 5     ,\
                                ybins    = 5     ,\
                                legend   = False ,\
                                rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "PSD2D"),\
                       self._extension)
            self._ph.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
