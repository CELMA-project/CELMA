#!/usr/bin/env python

"""Class for fourier modes plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap2
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotFourierModes
class PlotFourierModes(PlotSuperClass):
    """
    Class which contains the fourier modes data and the plotting configuration.
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

        # Set the member data
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, fourierModes, nModes = 7):
        #{{{docstring
        """
        Sets the fourier modes to be plotted.

        This function:
            * Sets the variable labels
            * Set the colors
            * Prepares the save name.

        Parameters
        ----------
        fourierModes : dict
            Dictionary where the keys are on the form "rho,z".
            The value is a dict containing of
            {varName:fourierModes, "time":time}
        nModes : int
            Number of modes to use
        """
        #}}}

        # Plus one as the offset mode will be excluded
        self._nModes  = nModes + 1

        # Set the member data
        self._fourierModes = fourierModes

        # Obtain the varname
        ind  = tuple(fourierModes.keys())[0]
        keys = fourierModes[ind].keys()
        self._varName = tuple(var for var in keys if var != "time")[0]
        # Strip the variable name if the
        self._varName.replace("Magnitude","")
        self._varName.replace("AngularVelocity","")

        # Obtain the color
        self._colors = seqCMap2(np.linspace(0, 1, len(fourierModes.keys())))

        self._prepareLabels()

        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)

        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}".format(self._varName, "fourierModes"))

        if self._extension is None:
            self._extension = "png"
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        #{{{docstring
        """
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
        self._varLabelTemplate = r"${{}}${}".format(unitsOrNormalization)

        # Set the time label
        self._timeLabel = self._ph.tTxtDict["tTxtLabel"]
    #}}}

    #{{{plotSaveShowFourierModes
    def plotSaveShowFourierModes(self):
        """
        Performs the actual plotting.
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sorted(self._fourierModes.keys())

        for key in keys:
            # Make the label
            rho, z = key.split(",")

            # Set values
            self._ph.rhoTxtDict  ["value"] =\
                    plotNumberFormatter(float(rho), None)
            self._ph.zTxtDict    ["value"] =\
                    plotNumberFormatter(float(z), None)

            rho   = self._ph.rhoTxtDict  ["constRhoTxt"].\
                            format(self._ph.rhoTxtDict)
            theta = self._ph.zTxtDict["constZTxt"].\
                            format(self._ph.zTxtDict),\

            # Make the const values
            title = (r"{}$,$ {}").format(rho, theta)

            # Range starts from one, as we excludes the offset mode
            for modeNr, color in zip(range(1, self._nModes), self._colors)
                label=r"$m_\theta={}$".format(modeNr)

                ax.plot(self._fourierModes[key]["time"],\
                        self._fourierModes[key][self._varName+"Magnitude"],\
                        color=color, label=label)

            # Set axis labels
            ax.set_xlabel(self._timeLabel)
            ax.set_ylabel(self._varLabel)

            # Set the title
            ax.set_title(title)

            # Make the plot look nice
            self._ph.makePlotPretty(ax, rotation = 45)

            if self._showPlot:
                plt.show()

            if self._savePlot:
                # Sets the save name
                self._fileName = "{}-rho{}-z{}.{}".\
                        format(self._fileName, rho, z, self._extension)
                self._ph.savePlot(fig, self._fileName)

            plt.close(fig)
    #}}}
#}}}
