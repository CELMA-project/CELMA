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

    #{{{__init___
    def __init__(self    ,\
                 *args   ,\
                 PDF     ,\
                 mode    ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets the member data
        3. Prepares the labels

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        PDF : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"pdfX":pdfX, "pdfY":"pdfY"}
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._PDF    = PDF
        self._colors = seqCMap3(np.linspace(0, 1, len(PDF.keys())))

        # Obtain the varname
        ind  = PDF.keys()[0]
        keys = PDF[ind].keys()
        self._varName = [var[:-4] for var in keys if "PDF" in var][0]

        # Set the labels
        # NOTE: The probability does not have any units, but in order to
        #       obtain the probability one has to integrate over
        #       the PDF.
        #       Thus will the PDF have the dimension of the inverse of
        #       what is on the x axis
        pltVarName   = self.ph.getVarPltName(self._varname)

        norm  = self.uc.conversionDict[self._varName]["normalization"]
        units = self.uc.conversionDict[self._varName]["units"]

        # Set the variable label
        if self.convertToPhysical:
            if mode == "normal":
                self._xLabel = r"${}$ $[{}]$"
            elif mode == "fluct":
                self._xLabel = r"$\tilde{{{}}}$ $[{}]$"
            self._xLabel = self._xLabel.format(pltVarName, units)
            self._yLabel = r"$\mathrm{{PDF}}(\tilde{{{}}})$".\
                    format(pltVarName)
        else:
            if mode == "normal":
                self._xLabel = r"${}{}$"
            elif mode == "fluct":
                self._xLabel = r"$\tilde{{{}}}{}$"
            self._xLabel = self._xLabel.format(pltVarName, norm)
            self._yLabel = r"$\mathrm{{PDF}}(\tilde{{{}}}{})$".\
                    format(pltVarName, norm)
    #}}}

    #{{{plotPDFs
    def plotPDFs(self):
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
