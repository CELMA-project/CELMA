#!/usr/bin/env python

"""Class for PDF plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotPDF
class PlotPDF(PlotSuperClass):
    """
    Class which contains the PDF data and the plotting configuration.
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
    def setData(self, PDF, mode):
        #{{{docstring
        """
        Sets the time traces to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        PDF : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"pdfX":pdfX, "pdfY":"pdfY"}
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        """
        #}}}

        # Set the member data
        self._PDF = PDF
        self._mode = mode

        # Obtain the varname
        ind  = list(PDF.keys())[0]
        keys = PDF[ind].keys()
        self._varName = [var[:-4] for var in keys if "PDF" in var][0]

        # Obtain the color (pad away brigthest colors)
        pad = 3
        self._colors = seqCMap3(np.linspace(0, 1, len(PDF.keys())+pad))

        self._prepareLabels()

        # Set the labels
        pltVarName = self._ph.getVarPltName(self._varName)

        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])
        self._xLabel = self._xLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, "PDF", self._fluctName))

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
        # NOTE: The probability does not have any units, but in order to
        #       obtain the probability one has to integrate over
        #       the PDF.
        #       Thus will the PDF have the dimension of the inverse of
        #       what is on the x axis
        pdfStr = r"\mathrm{{PDF}}("

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
        if self.uc.convertToPhysical:
            unitsOrNormalization       = "[1/{units}]"
            xLabelunitsOrNormalization = "[{units}]"
            self._varLabelTemplate = r"$"+\
                                     pdfStr +\
                                     normalOrFluct +\
                                     r")$" +\
                                     r" $" + unitsOrNormalization + r"$"
        else:
            unitsOrNormalization       = "{normalization}"
            xLabelunitsOrNormalization = "{normalization}"
            self._varLabelTemplate = r"$"+\
                                     pdfStr +\
                                     normalOrFluct +\
                                     unitsOrNormalization +\
                                     r")$"

        # Set the x-label
        self._xLabelTemplate = r"$" + normalOrFluct + r"$" +\
               r" $" + xLabelunitsOrNormalization + r"$"
    #}}}
    #{{{plotSaveShowPDF
    def plotSaveShowPDF(self):
        """ Plots the probability density function."""

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sorted(self._PDF.keys())

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

            ax.plot(self._PDF[key]["{}PDFX".format(self._varName)],\
                    self._PDF[key]["{}PDFY".format(self._varName)],\
                    color=color,\
                    label=label)

        # Use logarithmic scale
        ax.set_yscale("log")

        # Set axis labels
        ax.set_xlabel(self._xLabel)
        ax.set_ylabel(self._varLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(ax, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "PDF"),\
                       self._extension)
            self._ph.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
