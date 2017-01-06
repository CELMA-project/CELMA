#!/usr/bin/env python

"""Class for radialProfile plots"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import matplotlib.pylab as plt
import numpy as np
import os

# NOTE: May be suffering from a DRY case

#{{{PlotProfAndGradCompare
class PlotProfAndGradCompare(PlotSuperClass):
    """
    Class which contains the profile and gradient data and the plotting
    configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (10,16), **kwargs):
        #{{{docstring
        """
        Constructor for PlotProfAndGradCompare

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
    def setData(self, profGrad):
        #{{{docstring
        """
        Sets the profile and gradient to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        profGrad : dict
            Dictionary where the keys are:
                * varName      - Name of the variable under investigation
                * steadyVar    - The profile in the last steady state
                                 time point
                * varAvg       - The averaged profile
                * varAvgStd    - The standard deviation of the average
                * DDXSteadyVar - The radial derivative of the steady
                                 state path
                * DDXVarAvg    - The radial derivative of the average
                * DDXVarAvgStd - The standard deviation of the radial
                                 derivative of the average
                * X            - The radial coordinate
                * zPos         - The fixed z position
        """
        #}}}

        # Set the member data
        self._varName  = profGrad.pop("varName")
        self._rho      = profGrad.pop("X")
        self._zPos     = profGrad.pop("zPos")
        self._profGrad = profGrad

        # Obtain the color (pad away brigthest colors)
        pad = 1
        self._colors = seqCMap3(np.linspace(0, 1, 2+pad))
        self._lines  = ("-", "-.")

        # Prepare the labels
        self._prepareLabels()

        # Set the labels
        pltVarName    = self._ph.getVarPltName(self._varName)
        DDXPltVarName = r"\partial_\rho " + pltVarName

        self._profileLabel =\
                self._profileLabelTemplate.\
                format(pltVarName, **self.uc.conversionDict[self._varName])
        self._gradLabel =\
                self._gradLabelTemplate.\
                format(pltVarName, **self.uc.conversionDict[self._varName])
        self._gradAvgLabel =\
                self._avgLabelTemplate.\
                format(DDXPltVarName, **self.uc.conversionDict[self._varName])

        # Set the legends
        self._varAvgLegend =\
                self._avgLabelTemplate.\
                format(pltVarName, **self.uc.conversionDict[self._varName])
        self._varGradLegend =\
                self._varGradLegendTemplate.\
                format(pltVarName, **self.uc.conversionDict[self._varName])
        self._steadyStateLegend =\
                self._steadyStateLegendTemplate.format(pltVarName)
        self._steadyStateGradLegend =\
                self._steadyStateLegendTemplate.format(DDXPltVarName)

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}".format(self._varName, "profAndGradCompare"))

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

        self._profileLabelTemplate = r"${}$" + unitsOrNormalization
        self._gradLabelTemplate    = r"$\partial_\rho {}$" + gradUnitsOrNorm

        self._steadyStateLegendTemplate =\
            r"${}_{{\mathrm{{Steady \quad state}}}}$"
        self._avgLabelTemplate =\
            r"$\langle\langle {0}"+norm+r"\rangle_\theta\rangle_t$"
        self._varGradLegendTemplate =\
            r"$\partial_\rho $"+self._avgLabelTemplate

        # Set the title
        self._ph.zTxtDict    ["value"] =\
                plotNumberFormatter(float(self._zPos), None)
        self._title = "{}".format(self._ph.zTxtDict["constZTxt"].\
                                  format(self._ph.zTxtDict))

        # Set the x-axis label
        self._xLabel = self._ph.rhoTxtDict["rhoTxtLabel"]
    #}}}

    #{{{setDataForSecondVar
    def setDataForSecondVar(self, profGrad):
        #{{{docstring
        """
        Sets the additional fluctuation profile to be plotted.

        This function calls setData, and also alters the variable
        labels, colors and the save name.

        Parameters
        ----------
        profGrad : dict
            Dictionary where the keys are (in addition to those listed
            in setData:
                * varName2   - Name of the second variable under investigation
                * var2AvgStd - The standard deviation of the average of
                               the second variable
        """
        #}}}

        # Set the member data
        self._var2Name   = profGrad.pop("var2Name")
        self._var2AvgStd = profGrad["var2AvgStd"]

        # Obtain the color (pad away brigthest colors)
        pad = 1
        self._colors = seqCMap3(np.linspace(0, 1, 3+pad))
        self._lines  = ("--", "-", "-")

        # Prepare the labels
        self._prepareLabelsForSecondVar()

        # Set the labels
        pltVarName  = self._ph.getVarPltName(self._varName)
        pltVar2Name = self._ph.getVarPltName(self._var2Name)
        self._stdLabel1 =\
                self._stdLabelTemplate.\
                format(pltVarName, **self.uc.conversionDict[self._varName])
        self._stdLabel2 =\
                self._stdLabelTemplate.\
                format(pltVar2Name, **self.uc.conversionDict[self._var2Name])

        # Set the legends
        self._var1StdLegend =\
                self._stdLegendTemplate.\
                format(pltVarName, **self.uc.conversionDict[self._varName])
        self._var2StdLegend =\
                self._stdLegendTemplate.\
                format(pltVar2Name, **self.uc.conversionDict[self._var2Name])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, self._var2Name, "posOfFluct"))

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabelsForSecondVar
    def _prepareLabelsForSecondVar(self):
        """
        Prepares the labels for plotting.
        """

        # Set var label template
        if self.uc.convertToPhysical:
            norm = ""
            units = r" $[{units}]$"
        else:
            norm = r"${normalization}$"
            units = ""

        self._stdLegendTemplate = r"$\sqrt{{\langle\langle({0}"+norm+\
                                  r"-\langle\langle {0}"+norm+\
                                  r"\rangle_\theta\rangle_t)^2"+\
                                  r"\rangle_\theta\rangle_t}}$"

        self._stdLabelTemplate = r"${}$" + units + norm
    #}}}

    #{{{plotSaveShowRadialProfiles
    def plotSaveShowRadialProfiles(self):
        """
        Performs the actual plotting.

        Only setData needs to be called before calling this function
        """

        # Create the plot
        fig, (varAx, DDXVarAx) =\
                plt.subplots(nrows=2, figsize=self._pltSize, sharex=True)

        self._plotProfAndGrad(varAx, DDXVarAx)

        # Set the title
        varAx.set_title(self._title)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.1, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}

    #{{{_plotProfAndGrad
    def _plotProfAndGrad(self, varAx, DDXVarAx):
        #{{{docstring
        """
        Plots the profiles and the gradient

        Parameters
        ----------
        varAx : axis
            Axis to plot the profile on
        DDXVarAx : axis
            Axis to plot the gradient on
        """
        #}}}
        # Plot on the varAx
        # Steady state
        varAx.plot(self._rho,\
                   self._profGrad["steadyVar"],\
                   color     = self._colors[0],\
                   linestyle = self._lines[0],\
                   label     = self._steadyStateLegend
                  )
        # Average
        varAx.plot(self._rho,\
                   self._profGrad["varAvg"],\
                   color     = self._colors[1],\
                   linestyle = self._lines[1],\
                   label     = self._varAvgLegend
                  )
        # Fill
        varAx.fill_between(\
                self._rho,\
                self._profGrad["varAvg"]+self._profGrad["varAvgStd"],\
                self._profGrad["varAvg"]-self._profGrad["varAvgStd"],\
                facecolor=self._colors[1], edgecolor="none", alpha=0.2)

        # Set decorations
        varAx.set_ylabel(self._profileLabel)

        # Plot on the DDXAx
        # Steady state
        DDXVarAx.plot(self._rho,\
                      self._profGrad["DDXSteadyVar"],\
                      color     = self._colors[0],\
                      linestyle = self._lines[0],\
                      label     = self._steadyStateGradLegend
                     )
        # Average
        DDXVarAx.plot(self._rho,\
                      self._profGrad["DDXVarAvg"],\
                      color     = self._colors[1],\
                      linestyle = self._lines[1],\
                      label     = self._varGradLegend
                     )
        # Fill
        DDXVarAx.fill_between(\
                self._rho,\
                self._profGrad["DDXVarAvg"]+self._profGrad["DDXVarAvgStd"],\
                self._profGrad["DDXVarAvg"]-self._profGrad["DDXVarAvgStd"],\
                facecolor=self._colors[1], edgecolor="none", alpha=0.2)

        # Set decorations
        DDXVarAx.set_xlabel(self._xLabel)
        DDXVarAx.set_ylabel(self._gradLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(varAx, loc="upper right", ybins = 6)
        self._ph.makePlotPretty(DDXVarAx, loc="lower right",\
                                rotation = 45, ybins = 6)
        #}}}

    #{{{plotSaveShowPosOfFluct
    def plotSaveShowPosOfFluct(self):
        """
        Performs the actual plotting.

        setData and setDataForSecondVar needs to be called before
        calling this function.
        """

        # Create the plot
        fig, (varAx, DDXVarAx, fluct1Ax, fluct2Ax) =\
                plt.subplots(nrows=4, figsize=self._pltSize, sharex=True)

        self._plotProfAndGrad(varAx, DDXVarAx)

        # Plot the fluctuations
        fluct1Ax.plot(self._rho                      ,\
                      self._profGrad["varAvgStd"]    ,\
                      color     = self._colors[1]    ,\
                      linestyle = self._lines [1]    ,\
                      label     = self._var1StdLegend,\
                     )

        # Set decorations
        fluct1Ax.set_ylabel(self._stdLabel1)

        fluct2Ax.plot(self._rho                      ,\
                      self._var2AvgStd               ,\
                      color     = self._colors[2]    ,\
                      linestyle = self._lines [2]    ,\
                      label     = self._var2StdLegend,\
                     )

        # Set decorations
        fluct2Ax.set_ylabel(self._stdLabel2)
        DDXVarAx.set_xlabel("")
        fluct2Ax.set_xlabel(self._xLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(fluct1Ax, loc="lower right", ybins = 6)
        self._ph.makePlotPretty(fluct2Ax, loc="lower right",\
                                rotation = 45, ybins = 6)

        # Set the title
        varAx.set_title(self._title)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.2, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
