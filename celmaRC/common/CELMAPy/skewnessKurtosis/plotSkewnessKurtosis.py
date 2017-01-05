#!/usr/bin/env python

"""Class for skewness and kurtosis plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotSkewnessKurtosis
class PlotSkewnessKurtosis(PlotSuperClass):
    """
    Class which contains the skewness and kurtosis data and the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (9,5), **kwargs):
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
    def setData(self, skewKurt, mode):
        #{{{docstring
        """
        Sets the skewness and kurtosis to be plotted.

        This function also sets:
            * The variable labels
            * The colors
            * The save name

        Parameters
        ----------
        skewKurt : dict
            Dictionary with the keys
                * skew     - array of the skewnesses
                * kurt     - array of the kurtosises
                * rho      - array of the rho coordinate
                * thetaPos - integer with the fixed theta position
                * zPos     - integer with the fixed z position
                * varName  - the variabel name
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        """
        #}}}

        # Set the member data
        self._mode = mode
        self._skew = skewKurt.pop("skew")
        self._kurt = skewKurt.pop("kurt")
        self._rho  = skewKurt.pop("rho")

        # Obtain the varname
        self._varName = skewKurt.pop("varName")

        # Obtain the color (pad away brigthest colors)
        pad = 2
        self._colors = seqCMap3(np.linspace(0, 1, 2+pad))
        self._lines  = ("-", "-.")

        self._prepareLabels()

        # Set the labels
        pltVarName       = self._ph.getVarPltName(self._varName)
        self._skewLegend = self._skewLegendTemplate.format(pltVarName)
        self._kurtLegend = self._kurtLegendTemplate.format(pltVarName)

        # Set the title
        # Set values
        z     = skewKurt.pop("zPos")
        theta = skewKurt.pop("thetaPos")
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
                "{}-{}-{}".format(self._varName, "skewKurt", self._fluctName))

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
        # NOTE: Although the skewness and kurtosis has units, we chose
        #       here not to use them.

        skewStr = r"\mathrm{{Skewness}}("
        kurtStr = r"\mathrm{{Excess \quadkurtosis}}("

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

        self._skewLegendTemplate = r"$"+skewStr+normalOrFluct+r")$"
        self._kurtLegendTemplate = r"$"+kurtStr+normalOrFluct+r")$"

        # Set the x-axis label
        self._rhoLabel = self._ph.rhoTxtDict["rhoTxtLabel"]
    #}}}

    #{{{plotSaveShowSkewnessKurtosis
    def plotSaveShowSkewnessKurtosis(self):
        """
        Performs the actual plotting.
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        # Plot
        ax.plot(self._rho, self._skew,\
                linestyle = self._lines[0], color=self._colors[0],\
                label=self._skewLegend)
        ax.plot(self._rho, self._kurt,\
                linestyle = self._lines[1], color=self._colors[1],\
                label=self._kurtLegend)

        # Set axis label
        ax.set_xlabel(self._rhoLabel)

        # Set the title
        ax.set_title(self._title)

        # Make the plot look nice
        self._ph.makePlotPretty(ax, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath, "skewKurt"),\
                       self._extension)
            self._ph.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
