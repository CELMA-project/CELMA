#!/usr/bin/env python

"""Class for radial flux plot"""

from ..plotHelpers import plotNumberFormatter, seqCMap2, seqCMap3
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
import numpy as np
import os

YOU ARE HERE:
THESE ARE FIXMES
▸ commonDrivers/
▸ energy/
▸ fieldPlotters/
▸ growthRates/
▸ radialFlux/
▸ skewnessKurtosis/


#{{{PlotRadialFlux
class PlotRadialFlux(PlotsSuperClass):
    """
    Class which contains the radial flux data, and the plotting configuration.
    """

    #{{{__init___
    def __init__(self      ,\
                 *args     ,\
                 radialFlux,\
                 mode      ,\
                 **kwargs  ,\
                 ):
        #{{{docstring
        r"""
        This constructor:

        1. Calls the parent constructor
        2. Sets the member data
        3. Prepares the labels

        Parameters
        ----------
        radialFLux : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"pdfX":pdfX, "pdfY":"pdfY"}
        mode : ["normal"|"avg"|"fluct"]
            If mode is "normal" the output is on the form nu.
            If mode is "avg" the output is on the form <nu>.
            If mode is "fluct" the output is on the form \tilde{n}\tilde{u}.
            Note that
            <nu> = <(<n>+\tidle{n})(<u>+\tidle{u})>
                 = <<n><u>> + <\tidle{n}\tidle{u})>
                 = <n><u> + <\tidle{n}\tidle{u})>
            So that <n><u> is given by the "avg" - "fluct"
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._PDF    = PDF
        self._colors = seqCMap3(np.linspace(0, 1, len(timeTraces.keys())))

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
        if probes.helper.convertToPhysical:
            if mode == "normal":
                self._xLabel = r"${}$ $[{}]$"
            elif mode == "fluct"
                self._xLabel = r"$\tilde{{{}}}$ $[{}]$"
            self._xLabel = self._xLabel.format(pltVarName, Units)
            self._yLabel = r"$\mathrm{{PDF}}(\tilde{{{}}})$".\
                    format(pltVarName))
        else:
            if mode == "normal":
                self._xLabel = r"${}{}$"
            elif mode == "fluct"
                self._xLabel = r"$\tilde{{{}}}{}$"
            self._xLabel = self._xLabel.format(pltVarName, norm)
            self._yLabel = r"$\mathrm{{PDF}}(\tilde{{{}}}{})$".\
                    format(pltVarName, norm)
    #}}}

    #{{{plotPDFs
    def plotPDFs(self):
        """ Plots the time traces."""

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sort(self._timeTraces.keys())

        for key, color in keys, self._colors:
            # Make the label
            rhoInd, thetaInd, zInd = key.split(",")
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
