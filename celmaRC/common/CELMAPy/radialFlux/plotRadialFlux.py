#!/usr/bin/env python

"""Class for radial flux plot"""

from ..superClasses import PlotsSuperClass
from ..plotHelpers import seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

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
                 **kwargs):
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
            {"RadialFluxX":RadialFluxX, "RadialFluxY":"RadialFluxY"}
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
        self._radialFlux    = radialFlux
        self._colors = seqCMap3(np.linspace(0, 1, len(radialFlux.keys())))

        # Obtain the varname
        ind  = radialFlux.keys()[0]
        keys = radialFlux[ind].keys()
        self._varName = [var[:-4] for var in keys if "RadialFlux" in var][0]

        # Set the labels
        # NOTE: The probability does not have any units, but in order to
        #       obtain the probability one has to integrate over
        #       the RadialFlux.
        #       Thus will the RadialFlux have the dimension of the inverse of
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
            self._yLabel = r"$\mathrm{{RadialFlux}}(\tilde{{{}}})$".\
                    format(pltVarName)
        else:
            if mode == "normal":
                self._xLabel = r"${}{}$"
            elif mode == "fluct":
                self._xLabel = r"$\tilde{{{}}}{}$"
            self._xLabel = self._xLabel.format(pltVarName, norm)
            self._yLabel = r"$\mathrm{{RadialFlux}}(\tilde{{{}}}{})$".\
                    format(pltVarName, norm)
    #}}}

    #{{{plotRadialFluxes
    def plotRadialFluxes(self):
        """ Plots the radial fluxes."""

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sorted(self._timeTraces.keys())

        for key, color in keys, self._colors:
            # Make the label
            rho, theta, z = key.split(",")
            label = (r"$\rho={}$ $\theta={}$ $z={}$").format(rho, theta. z)

            ax.plot(self._RadialFlux[key]["time"],\
                    self._RadialFlux[key][self._varName],\
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
                format(os.path.join(self._savePath, "radialFlux"),\
                       self._extension)
            self.ph.savePlot(fig, fileName, (self._leg,))

        plt.close(fig)
    #}}}
#}}}
