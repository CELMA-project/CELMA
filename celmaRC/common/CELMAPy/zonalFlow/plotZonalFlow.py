#!/usr/bin/env python

"""
Class for zonalFlow plots
"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import matplotlib.pylab as plt
import numpy as np
import os

# NOTE: May be suffering from a DRY case

#{{{PlotzonalFlow
class PlotzonalFlow(PlotSuperClass):
    """
    Class which contains the profile and gradient data and the plotting
    configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (10,16), **kwargs):
        #{{{docstring
        """
        Constructor for PlotzonalFlow

        * Calls the parent class
        * Sets the spatial title
        * Sets the xlabel

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        pltSize : tuple
            The size of the plot
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, polExB):
        #{{{docstring
        """
        Sets the profiles and gradients to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        polExB : dict
            Dictionary where the keys are:
                * polExBSS          - The poloidal ExB velocity profile
                                      in the last steady state time point
                * angPolExBSS       - The angular poloidal ExB velocity profile
                                      in the last steady state time point
                * angPolExBShearSS  - The shear of the angular poloidal ExB
                                      velocity profile in the steady
                                      state time point
                * polExBAvg         - The average poloidal ExB velocity
                                      profile in the turbulent state
                * polExBStd         - The standard deviation of the poloidal
                                      ExB velocity profile in the
                                      turbulent state
                * angPolExBAvg      - The average angluar poloidal ExB velocity
                                      profile in the turbulent state
                * angPolExBStd      - The standard deviation of the
                                      angluar poloidal ExB velocity
                                      profile in the turbulent state
                * angPolExBShearAvg - The average shear of the angular poloidal
                                      ExB velocity profile int the
                                      turbulent state
                * angPolExBShearStd - The standard deviation of the shear of
                                      the angular poloidal ExB velocity
                                      profile int the turbulent state
                * rho               - The radial coordinate
                * zPos              - The fixed z position
        """
        #}}}

        # Set the member data
        self._rho    = polExB.pop("rho")
        self._zPos   = polExB.pop("zPos")
        self._polExB = polExB

        # Obtain the color (pad away brigthest colors)
        pad = 1
        self._colors = seqCMap3(np.linspace(0, 1, 2+pad))
        self._lines  = ("-", "-.")

        # Prepare the labels
        self._prepareLabels()

        # Set the labels
        # NOTE: We use growth rate to find the frequency
        self._polExBLabel = self._polExBLabelTemplate.\
                format(**self.uc.conversionDict["u"])
        self._angPolExBLabel = self._angPolExBLabelTemplate.\
                format(**self.uc.conversionDict["growthRate"])
        self._angGradPolExBLabel = self._angGradPolExBLabelTemplate.\
                format(**self.uc.conversionDict["growthRate"])

        # Set the legends
        # NOTE: We use format to get away the double {{
        # Seady state
        self._sSPolExBLegend       =self._sSPolExBLegend       .format()
        self._sSAngPolExBLegend    =self._sSAngPolExBLegend    .format()
        self._gradSSAngPolExBLegend=self._gradSSAngPolExBLegend.format()

        # Turbulent
        self._polExBLegend = self._polExBLegendTemplate.\
                format(**self.uc.conversionDict["u"])
        self._angPolExBLegend = self._angPolExBLegendTemplate.\
                format(**self.uc.conversionDict["growthRate"])
        self._gradAngPolExBLegend = self._gradAngPolExBLegendTemplate.\
                format(**self.uc.conversionDict["growthRate"])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}".format(self._varName, "zonalFlow"))

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Set var label template
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
            norm                 = ""
            gradUnitsOrNorm      = r" $[{units}\mathrm{{m}}^{{-1}}]$"
        else:

            unitsOrNormalization = r"${normalization}$"
            norm                 = r"${normalization}$"
            # NOTE: The normalization always have a divided by before
            #       the end of the string. Therefore the added rho_s
            #       will appear as a "divided by"
            gradUnitsOrNorm      = r"${normalization}/\rho_s$"

        # Set the labels
        polExBStr    = r"$u_{{\mathrm{{E\times B}}, \theta}}$"
        angPolExBStr = r"$\omega_{{\mathrm{{E\times B}}, \theta}}$"
        gradStr      = r"$\partial_\rho$"
        self._polExBLabelTemplate        = polExBStr + unitsOrNormalization
        self._angPolExBLabelTemplate     = angPolExBStr + unitsOrNormalization
        self._angGradPolExBLabelTemplate = gradStr      +\
                                           angPolExBStr +\
                                           unitsOrNormalization

        # Set the legends
        # Seady state
        sSStr = r"$_{{\mathrm{{Steady \quad state}}}}$"
        self._sSPolExBLegend        = polExBStr + sSStr
        self._sSAngPolExBLegend     = angPolExBStr + sSStr
        self._gradSSAngPolExBLegend = grad + StrangPolExBStr + sSStr

        # Turbulent
        lAngle = r"$\langle\langle$"
        rAngle = r"$\rangle_\theta\rangle_t$"
        self._polExBLegendTemplate        =\
            lAngle + self._sSPolExBLegend + norm + rAngle
        self._angPolExBLegendTemplate     =\
            lAngle + self._gradSSAngPolExBLegend + norm + rAngle
        self._gradAngPolExBLegendTemplate =\
            lAngle + self._sSAngPolExBLegend + norm+ r Angle

        # Set the title
        self._ph.zTxtDict    ["value"] =\
                plotNumberFormatter(float(self._zPos), None)
        self._title = "{}".format(self._ph.zTxtDict["constZTxt"].\
                                  format(self._ph.zTxtDict))

        # Set the x-axis label
        self._xLabel = self._ph.rhoTxtDict["rhoTxtLabel"]
    #}}}

    #{{{plotSaveShowPosOfFluct
    def plotSaveShowPosOfFluct(self):
        """
        Performs the actual plotting.

        setData and setDataForSecondVar needs to be called before
        calling this function.
        """

        # Create the plot
        fig, (polExBAx, angPolExBAx, shearPolExBAx) =\
                plt.subplots(nrows=3, figsize=self._pltSize, sharex=True)

        # Plot on the polExBAx
        # Steady state
        polExBAx.plot(self._rho,\
                       self._polExB["polExBSS"],\
                       color     = self._colors[0],\
                       linestyle = self._lines[0],\
                       label     = self._sSPolExBLegend
                      )
        # Average
        polExBAx.plot(self._rho,\
                      self._polExB["polExBAvg"],\
                      color     = self._colors[1],\
                      linestyle = self._lines[1],\
                      label     = self._polExBLegend
                     )
        # Fill
        polExBAx.fill_between(\
                    self._rho,\
                    self._polExB["polExBAvg"]+self._polExB["polExBStd"],\
                    self._polExB["polExBAvg"]-self._polExB["polExBStd"],\
                    facecolor=self._colors[1], edgecolor="none", alpha=0.2)

        # Set decorations
        polExBAx.set_ylabel(self._polExBLabel)

        # Plot on the angPolExBAx
        # Steady state
        angPolExBAx.plot(self._rho,\
                         self._polExB["angPolExBSS"],\
                         color     = self._colors[0],\
                         linestyle = self._lines[0],\
                         label     = self._sSAngPolExBLegend
                        )
        # Average
        angPolExBAx.plot(self._rho,\
                         self._polExB["angPolExBAvg"],\
                         color     = self._colors[1],\
                         linestyle = self._lines[1],\
                         label     = self._angPolExBLegend
                        )
        # Fill
        angPolExBAx.fill_between(\
                    self._rho,\
                    self._polExB["angPolExBAvg"]+self._polExB["angPolExBStd"],\
                    self._polExB["angPolExBAvg"]-self._polExB["angPolExBStd"],\
                    facecolor=self._colors[1], edgecolor="none", alpha=0.2)

        # Set decorations
        angPolExBAx.set_ylabel(self._angPolExBLabel)

        # Plot on the shearAngPolExBAx
        # Steady state
        shearAngPolExBAx.plot(self._rho,\
                               self._polExB["angPolExBShearSS"],\
                               color     = self._colors[0],\
                               linestyle = self._lines[0],\
                               label     = self._gradSSAngPolExBLegend
                              )
        # Average
        angPolExBAx.plot(self._rho,\
                         self._polExB["angPolExBShearAvg"],\
                         color     = self._colors[1],\
                         linestyle = self._lines[1],\
                         label     = self._gradAngPolExBLegend
                        )
        # Fill
        angPolExBAx.fill_between(\
           self._rho,\
           self._polExB["angPolExBShearAvg"]+self._polExB["angPolExBshearStd"],\
           self._polExB["angPolExBShearAvg"]-self._polExB["angPolExBshearStd"],\
           facecolor=self._colors[1], edgecolor="none", alpha=0.2)

        # Set decorations
        angPolExBAx.set_ylabel(self._angGradPolExBLabel)
        angPolExBAx.set_xlabel(self._xLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(polExBAx        , loc="lower right", ybins = 6)
        self._ph.makePlotPretty(angPolExBAx     , loc="lower right", ybins = 6)
        self._ph.makePlotPretty(shearAngPolExBAx, loc="lower right",\
                                rotation = 45, ybins = 6)

        # Set the title
        polExBAx.set_title(self._title)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.2, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
