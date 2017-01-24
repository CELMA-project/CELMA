#!/usr/bin/env python

"""Class for magnitude spectrum plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap2
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotMagnitudeSpectrum
class PlotMagnitudeSpectrum(PlotSuperClass):
    """
    Class which contains the magnitude spectrum data and the plotting configuration.
    """
    #{{{Static members
    _errorbarOptions = {"color"     :"k",\
                        "fmt"       :"o",\
                        "fmt"       :"o",\
                        "markersize":10 ,\
                        "ecolor"    :"k",\
                        "capsize"   :7  ,\
                        "elinewidth":3  ,\
                        }

    _markeredgewidth = 3
    #}}}

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

        # Set the member data
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, magnitudeSpectrum, varName):
        #{{{docstring
        """
        Sets the magnitude spectrum to be plotted.

        This function:
            * Sets the variable labels
            * Set the colors
            * Prepares the save name.

        Parameters
        ----------
        magnitudeSpectrum : dict
            The dictionary containing the mode spectra.
            The keys of the dict is on the form (rho,z), with a new dict
            as the value. This dict contains the following keys:
                * "modeAvg" - The average magnitude of the mode
                * "modeStd" - The standard deviation of the magnitude of
                              the mode
                * "modeNr"  - The mode number
        varName : str
            Name of the variable.
        """
        #}}}

        # Set the member data
        self._mSpec   = magnitudeSpectrum
        self._varName = varName

        self._prepareLabels()

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}".format(self._varName, "magnitudeSpectrum"))

        if self._extension is None:
            self._extension = "png"
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        #{{{docstring
        r"""
        Prepares the labels for plotting.

        NOTE: As we are taking the fourier transform in the
              dimensionless theta direction, that means that the units
              will remain the same as

              \hat{f}(x) = \int_\infty^\infty f(x) \exp(-i2\pi x\xi) d\xi
        """
        #}}}

        # Set var label template
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
        else:
            unitsOrNormalization = "${normalization}$"

        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)

        self._varLabelTemplate = r"${{}}${}".format(unitsOrNormalization)

        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])

        # FIXME: This is just a dirty way of doing things, consider
        #        refactor for clearity
        if self.uc.convertToPhysical:
            theSplit = self._varLabel.split("[")
            # Exclude the last $
            var = theSplit[0][:-1].replace("$", "")
            units = ("[" + theSplit[1]).replace("$", "")
            self._varLabel = r"$|\mathrm{{FT}}({})|$ ${}$".format(var, units)
        else:
            self._varLabel = r"$|\mathrm{{FT}}({})|$".format(self._varLabel)

        # Set the x-label
        self._xLabel = r"$\mathrm{Mode\;number}$"
    #}}}

    #{{{plotSaveShowMagnitudeSpectrum
    def plotSaveShowMagnitudeSpectrum(self):
        """
        Performs the actual plotting.
        """

        keys = sorted(self._mSpec.keys())

        for key in keys:
            # Create the plot
            fig = plt.figure(figsize = self._pltSize)
            ax  = fig.add_subplot(111)

            # Make the label
            rho, z = key.split(",")

            # Set values
            self._ph.rhoTxtDict  ["value"] =\
                    plotNumberFormatter(float(rho), None)
            self._ph.zTxtDict    ["value"] =\
                    plotNumberFormatter(float(z), None)

            rhoVal   = self._ph.rhoTxtDict["value"].replace("$","")
            rhoTitle = self._ph.rhoTxtDict  ["constRhoTxt"].\
                            format(self._ph.rhoTxtDict)
            zVal     = self._ph.zTxtDict["value"].replace("$","")
            zTitle   = self._ph.zTxtDict["constZTxt"].\
                            format(self._ph.zTxtDict)

            # Make the const values
            title = (r"{}$,$ {}").format(rhoTitle, zTitle)

            # Range starts from one, as we excludes the offset mode
            (_, caps, _) = ax.errorbar(\
                self._mSpec[key]["modeNr"]          ,\
                self._mSpec[key]["modeAvg"]         ,\
                yerr  = self._mSpec[key]["modeStd"] ,\
                **self._errorbarOptions)

            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            # Set axis labels
            ax.set_xlabel(self._xLabel)
            ax.set_ylabel(self._varLabel)

            # Set logarithmic scale
            ax.set_yscale("log")

            # Set the title
            ax.set_title(title)

            # Make the plot look nice
            self._ph.makePlotPretty(ax, rotation = 45, legend = False)

            if self._showPlot:
                plt.show()

            if self._savePlot:
                # Sets the save name
                fileName = "{}-rho-{}-z-{}.{}".\
                    format(self._fileName, rhoVal, zVal, self._extension)
                self._ph.savePlot(fig, fileName)

            plt.close(fig)
    #}}}
#}}}
