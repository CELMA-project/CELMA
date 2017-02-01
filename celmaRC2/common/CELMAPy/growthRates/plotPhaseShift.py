#!/usr/bin/env python

"""
Contains class to plot the phase shift
"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import SizeMaker, plotNumberFormatter, seqCMap2
import numpy as np
import matplotlib.pylab as plt
import os

#{{{PlotPhaseShift
class PlotPhaseShift(PlotSuperClass):
    """Class which contains the phase shift and the plotting configuration."""

    #{{{Static members
    _errorbarOptions = {"color"     :"k",\
                        "fmt"       :"o",\
                        "markersize":7  ,\
                        "ecolor"    :"k",\
                        "capsize"   :7  ,\
                        "capthick"  :3  ,\
                        "elinewidth":3  ,\
                        "alpha"     :0.7,\
                        }
    #}}}

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)
    #}}}

    #{{{setData
    def setData(self,\
                analyticPhaseShiftDF,\
                simulationPhaseShiftDF,\
                positionTuple):
        #{{{docstring
        """
        Sets the data to be plotted.

        This function:
            * Sets the data frame
            * Sets the variable labels
            * Sets the colors
            * Prepares the save name

        Parameters
        ----------
        analyticPhaseShiftDF : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "phaseShift" - Phase shift obtained from the
                                 analytical expression
            over the observation "modeNr" over the observation "Scan"
        simulationPhaseShiftDF : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "phaseShift" - Phase shift obtained from the simulations.
            over the observation "Scan"
        positionTuple : tuple
            The tuple containing (rho, theta, z).
        """
        #}}}

        # Set the data frames
        self._aPDF = analyticPhaseShiftDF
        self._sPDF = simulationPhaseShiftDF

        # Set the position
        rho, theta, z = positionTuple
        # Set values
        self._ph.rhoTxtDict  ["value"] = plotNumberFormatter(float(rho)  , None)
        self._ph.thetaTxtDict["value"] = plotNumberFormatter(float(theta), None)
        self._ph.zTxtDict    ["value"] = plotNumberFormatter(float(z)    , None)
        # Get the values
        rhoVal   = self._ph.rhoTxtDict  ["value"].replace("$","")
        thetaVal = self._ph.thetaTxtDict["value"].replace("$","")
        zVal     = self._ph.zTxtDict    ["value"].replace("$","")
        # Get the titles
        rhoTitle = self._ph.rhoTxtDict["constRhoTxt"].\
                    format(self._ph.rhoTxtDict)
        thetaTitle = self._ph.thetaTxtDict["constThetaTxt"].\
                    format(self._ph.thetaTxtDict)
        zTitle     = self._ph.zTxtDict  ["constZTxt"].\
                    format(self._ph.zTxtDict)
        self._title = r"{}$,$ {}$,$ {}".format(rhoTitle, thetaTitle, zTitle)


        self._prepareLabels()

        # Set colors
        self._colors =\
            seqCMap2(np.linspace(0, 1, len(self._aPDF.index.levels[0].values)))

        # Set the fileName
        fileName = "phaseShift"
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-rho-{}-theta-{}-z-{}".\
                format(fileName, rhoVal, thetaVal, zVal))

        if self._extension is None:
            self._extension = "png"
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        self._yLabel = r"$\Psi$"

        # Set var label templates
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
        else:
            unitsOrNormalization = "${normalization}$"
        self._varLabelTemplate = r"{{}}{}".format(unitsOrNormalization)

        # Get the plot decoration for the scan
        scanName    = self._aPDF.index.names[-1]
        scanPltName = self._ph.getVarPltName(scanName)

        # Get the x label
        self._xLabel = self._varLabelTemplate.\
                format(scanPltName, **self.uc.conversionDict[scanName])
    #}}}

    #{{{plotSaveShowPhaseShift
    def plotSaveShowPhaseShift(self):
        """
        Performs the plotting by calling the workers.
        """

        # Create the axes
        figSize = SizeMaker.standard(s=0.7)
        fig, ax = plt.subplots(figsize=figSize)

        # Plot the analytical phase shift
        for outerInd, color in\
                zip(self._aPDF.index.levels[0].values, self._colors):
            # Update the colors
            self._errorbarOptions.update({"color":color})

            ax.errorbar(\
                self._aPDF.loc[outerInd]["phaseShift"].index.values,\
                self._aPDF.loc[outerInd]["phaseShift"].values      ,\
                label = r"$m_\theta={}$".format(outerInd),\
                **self._errorbarOptions)

        # Plot the phase shift from the simulations
        # Update the colors and markers
        self._errorbarOptions.update({"color":"k", "fmt":"v"})

        ax.errorbar(\
            self._sPDF["phaseShift"].index.values,\
            self._sPDF["phaseShift"].values      ,\
            label = "$\mathrm{Simulation}$"      ,\
            **self._errorbarOptions)

        # Set labels
        ax.set_ylabel(self._yLabel)
        ax.set_xlabel(self._xLabel)

        # Tweak the ticks
        # Add 10% margins for readability
        ax.margins(x=0.1, y=0.1)

        # Make the plot look nice
        self._ph.makePlotPretty(ax, loc ="upper left",\
                                legend = True, rotation = 45)

        # Set the ticks
        ticks = tuple(self._aPDF.loc[outerInd]["phaseShift"].index.values)
        ax.xaxis.set_ticks(ticks)
        self._ph.radiansOnAxis(ax, "y")


        # Set the title
        fig.suptitle(self._title, y=0.9, transform=ax.transAxes, va="bottom")

        # Remove old legend
        leg = ax.legend()
        leg.remove()
        # Put the legend outside
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles,\
                   labels,\
                   bbox_to_anchor=(1.05, 1.0),\
                   loc="upper left",\
                   borderaxespad=0.,\
                   bbox_transform = ax.transAxes,\
                   )

        if self._showPlot:
            plt.show()

        if self._savePlot:
            # Sets the save name
            fileName = "{}-{}.{}".\
                format(self._fileName, self._aPDF.index.names[0], self._extension)

            self._ph.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
