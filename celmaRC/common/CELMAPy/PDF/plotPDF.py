#!/usr/bin/env python

"""Class for PDF plot"""

from ..superClasses import PlotsSuperClass
from ..plotHelpers import seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotPDF
class PlotPDF(PlotsSuperClass):
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
        ind  = PDF.keys()[0]
        keys = PDF[ind].keys()
        self._varName = [var[:-4] for var in keys if "PDF" in var][0]

        # Obtain the color (pad away brigthest colors)
        pad = 3
        self._colors = seqCMap3(np.linspace(0, 1, len(PDF.keys())+pad))

        self._prepareLabels()

        # Set the labels
        pltVarName = self._ph.getVarPltName(self._varName)

        # Make a unitsOrNormDict which can be altered
        # If not(self.uc.convertToPhysical) we need to manually modify
        # the normalization argument
        unitsOrNormDict = self.uc.conversionDict[self._varName]
        if not(self.uc.convertToPhysical):
            if r"/" in unitsOrNormDict["normalization"]:
                if unitsOrNormDict["normalization"][0] == r"/":
                    unitsOrNormDict["normalization"] =\
                        unitsOrNormDict["normalization"][1:]
                else:
                    # Swap the nominator and denominator
                    tmp = unitsOrNormDict["normalization"].split(r"/")
                    unitsOrNormDict["normalization"] =\
                            r"/".join((tmp[1],tmp[0]))
            else:
                # Add division sign
                unitsOrNormDict["normalization"] =\
                        r"/"+unitsOrNormDict["normalization"]

        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **unitsOrNormDict)

        self._xLabel = self._varLabelTemplate.\
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
        if self.uc.convertToPhysical:
            unitsOrNormalization       = " $[{units}]$"
            xLabelunitsOrNormalization = " $[1/{units}]$"
        else:
            unitsOrNormalization       = "${normalization}$"
            # NOTE: this will be treated at a later stage
            xLabelunitsOrNormalization = "${normalization}$"
        if self._mode == "normal":
            self._varLabelTemplate =\
                r"$\mathrm{{{PDF}}}({{}})${}".format(unitsOrNormalization)
            self._xLabel = r"${{}}${}".format(unitsOrNormalization)
            self._fluctName = ""
        elif self._mode == "fluct":
            self._varLabelTemplate =\
                r"$\mathrm{{{{{{{PDF}}}}}}}(\widetilde{{{{{{}}}}}})${}".\
                format(unitsOrNormalization)
            self._xLabel = r"$\widetilde{{{{{{}}}}}}${}".\
                               format(xLabelunitsOrNormalization)
            self._fluctName = "fluct"
        else:
            raise NotImplementedError("'{}'-mode not implemented.")
    #}}}

    #{{{plotSaveShowPDF
    def plotSaveShowPDF(self):
        """ Plots the probability density function."""

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sorted(self._PDF.keys())

        for key, color in keys, self._colors:
            # Make the label
            rho, theta, z = key.split(",")
            label = (r"$\rho={}$ $\theta={}$ $z={}$").format(rho, theta. z)

            ax.plot(self._PDF[key]["{}PDFX"].self._varName,\
                    self._PDF[key]["{}PDFY"].self._varName,\
                    color=color,\
                    label=label)

        # Set axis labels
        ax.set_xlabel(self._xLabel)
        ax.set_ylabel(self._yLabel)

        # Make the plot look nice
        self.ph.makePlotPretty(ax, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "PDF"),\
                       self._extension)
            self.ph.savePlot(fig, fileName, (self._leg,))

        plt.close(fig)
    #}}}
#}}}
