#!/usr/bin/env python

"""Class for PSD plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotPSD
class PlotPSD(PlotSuperClass):
    """
    Class which contains the PSD data and the plotting configuration.
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
    def setData(self, PSD, mode):
        #{{{docstring
        """
        Sets the time traces to be plotted.

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
            self._xLabel         = r"$1/t$ $[Hz]$"
        else:
            self._varLabelTemplate = r"$"+\
                                     psdStr +\
                                     normalOrFluct +\
                                     "{normalization}" +\
                                     r")$" +\
                                     r" $" + units + r"$"
            self._xLabel         = r"$1/t \omega_{ci}$"
    #}}}

    #{{{plotSaveShowPSD
    def plotSaveShowPSD(self):
        """
        Plots the probability density function.

        NOTE: The desibel use here is dividing by its own max, but not
              the global max.
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
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
            yVals = np.log10(       self._PSD[key]["{}PSDY".format(self._varName)][1:]/\
                             np.max(self._PSD[key]["{}PSDY".format(self._varName)][1:]))

            #Plot
            ax.plot(xVals, yVals, color=color, label=label)

        # Use logarithmic scale
        ax.set_xscale("log")

        # Set axis labels
        ax.set_xlabel(self._xLabel)
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
#}}}
