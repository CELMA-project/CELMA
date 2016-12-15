#!/usr/bin/env python

"""Class for PSD plot"""

from ..superClasses import PlotsSuperClass
from ..plotHelpers import seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotPSD
class PlotPSD(PlotsSuperClass):
    """
    Class which contains the PSD data and the plotting configuration.
    """

    #{{{__init___
    def __init__(self    ,\
                 *args   ,\
                 PSD     ,\
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
        PSD : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varPSDX:psdX, varPSDY:psdY}
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._PSD    = PSD
        self._colors = seqCMap3(np.linspace(0, 1, len(PSD.keys())))

        # Obtain the varname
        ind  = PSD.keys()[0]
        keys = PSD[ind].keys()
        self._varName = [var[:-4] for var in keys if "PSD" in var][0]

        # Set the labels
        # NOTE: The units will be in variableUnits**2/Hz for
        #       non-normalized variables.
        #       The normalization would be
        #       variableNormalization**2/(tOmegaCI)
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
            self._yLabel = r"$\mathrm{{PSD}}(\tilde{{{}}})$".\
                    format(pltVarName)
        else:
            if mode == "normal":
                self._xLabel = r"${}{}$"
            elif mode == "fluct":
                self._xLabel = r"$\tilde{{{}}}{}$"
            self._xLabel = self._xLabel.format(pltVarName, norm)
            self._yLabel = r"$\mathrm{{PSD}}(\tilde{{{}}}{})$".\
                    format(pltVarName, norm)
    #}}}

    #{{{plotPSDs
    def plotPSDs(self, mode="normal"):
        """
        Plots the time power spectral density.

        Parameters
        -----------
        mode : ["normal"|"dB"]
            The y-axis will be on logarithmic scale when mode is
            "normal", and on a log-log scale when the mode is "dB"
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sorted(self._PSD.keys())

        if mode == "normal":
            for key, color in keys, self._colors:
                # Make the label
                rho, theta, z = key.split(",")
                label = (r"$\rho={}$ $\theta={}$ $z={}$").format(rho, theta. z)
                
                # Clip the very first point as this is rediculously low
                ax.plot(self._probes.results[key]["psdX"][1:],\
                        self._probes.results[key]["psdY"][1:],\
                        color=color,\
                        label=label,\
                        alpha=self._alpha)

        elif mode == "dB":
            # Find max:
            curMax = 0
            for key, color in keys, self._colors:
                if np.max(self._probes.results[key]["psdY"][1:]) > curMax:
                    curMax = np.max(self._probes.results[key]["psdY"][1:])

            # Make the plots
            for key, color in keys, self._colors:
                # Make the label
                rho, theta, z = key.split(",")
                label = (r"$\rho={}$ $\theta={}$ $z={}$").format(rho, theta. z)
                
                # Clip the very first point as this is rediculously low
                ax.plot(self._probes.results[key]["psdX"][1:],\
                        np.log10(\
                            self._probes.results[key]["psdY"][1:]/\
                            curMax),\
                        color=color,\
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
    #}}}
#}}}
