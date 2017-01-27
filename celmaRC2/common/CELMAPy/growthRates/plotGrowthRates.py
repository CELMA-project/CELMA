#!/usr/bin/env python

"""
Contains class to plot the growth rates
"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import SizeMaker, plotNumberFormatter, seqCMap2, seqCMap3
import numpy as np
import matplotlib.pylab as plt
import os

#{{{PlotGrowthRates
class PlotGrowthRates(PlotSuperClass):
    """Class which contains the growth rates and the plotting configuration."""

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
                varName,\
                growthRateDataFrame,\
                positionTuple,\
                analytic = False):
        #{{{docstring
        """
        Sets the data to be plotted.

        This function:
            * Sets the data frame
            * Sets the variable labels
            * Prepares the save name

        Parameters
        ----------
        varName : str
            Name of the variable
        growthRateDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "growthRate"                  - Always present
                * "growthRateStd"               - Only if analytic = False
                * "averageAngularVelocity"      - Only if analytic = False
                * "averageAngularVelocityStd"   - Only if analytic = False
                * "angularVelocity"             - Only if analytic = True
            over the observation "modeNr" over the observation "Scan"
        positionTuple : tuple
            The tuple containing (rho, z).
        analytic : bool
            Wheter or not an analytical or semi-analytical expression
            has been used.
        """
        #}}}

        self._analytic = analytic

        # NOTE: Theta remains unspecified as we have done a fourier transform
        rho, z = positionTuple
        # Set values
        self._ph.rhoTxtDict["value"] = plotNumberFormatter(float(rho), None)
        self._ph.zTxtDict  ["value"] = plotNumberFormatter(float(z)  , None)
        # Get the values
        rhoVal= self._ph.rhoTxtDict["value"].replace("$","")
        zVal  = self._ph.zTxtDict  ["value"].replace("$","")
        # Get the titles
        rhoTitle=self._ph.rhoTxtDict["constRhoTxt"].format(self._ph.rhoTxtDict)
        zTitle  =self._ph.zTxtDict  ["constZTxt"]  .format(self._ph.zTxtDict)
        self._title = r"{}$,$ {}".format(rhoTitle, zTitle)

        # Set the growth rate
        self._gRDF = growthRateDataFrame

        # The variable name is stored in the outermost level of the data
        # frame observations
        self._varName = varName

        if not(analytic):
            # Set the plot name
            self._pltVarName = self._ph.getVarPltName(self._varName)

        self._prepareLabels()

        # Set the fileName
        grName = "growthRatesAnalytic" if analytic else "growthAnalytic"
        fileVarName = self._varName if analytic else self._varName + "-"
        self._fileName =\
            os.path.join(self._savePath,\
                "{}{}-rho-{}-z-{}".\
                format(fileVarName, grName, rhoVal, zVal))

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

        # Set var label and legend templates
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
            self._legendTemplate = r"{0[scanOrMode]} = {2[val]}${1[units]}$"
        else:
            unitsOrNormalization = "${normalization}$"
            self._legendTemplate =\
                r"{0[scanOrMode]}{1[normalization]} = ${2[val]}$"
        self._varLabelTemplate = r"{{}}{}".format(unitsOrNormalization)

        # Set the y-axes
        if not(self._analytic):
            imag = r"$\Im(\omega_{})$".format(self._pltVarName.replace("$",""))
            real = r"$\Re(\omega_{})$".format(self._pltVarName.replace("$",""))
        else:
            imag = r"$\Im(\omega)$"
            real = r"$\Re(\omega)$"
        self._imLabel =\
            self._varLabelTemplate.\
                format(imag, **self.uc.conversionDict["growthRate"])
        self._reLabel =\
            self._varLabelTemplate.\
                format(real, **self.uc.conversionDict["growthRate"])

        # Get the plot decoration for the scan
        self._scan         = self._gRDF.index.names[0]
        scanName           = self._ph.getVarPltName(self._scan)
        self._scanNameDict = {"scanOrMode":scanName}
        self._scanNormalizationOrUnitsDict = self.uc.conversionDict[self._scan]

        # Get the plot decoration for the modes
        self._mode         = self._gRDF.index.names[1]
        modeName           = self._ph.getVarPltName(self._mode)
        self._modeNameDict = {"scanOrMode":modeName}
        self._modeNormalizationOrUnitsDict = self.uc.conversionDict[self._mode]

        # Get the x labels
        self._scanXLabel = self._varLabelTemplate.\
                format(modeName, **self.uc.conversionDict[self._mode])
        self._modeXLabel = self._varLabelTemplate.\
                format(scanName, **self.uc.conversionDict[self._scan])
    #}}}

    #{{{plotSaveShowGrowthRates
    def plotSaveShowGrowthRates(self):
        """
        Performs the plotting by calling the workers.
        """

        if not(self._analytic):
            self._gRDF.rename(columns=\
                    {"averageAngularVelocity":"angularVelocity",\
                     "averageAngularVelocityStd":"angularVelocityStd"})

        # Make the scan plot
        scanColors =\
            seqCMap3(np.linspace(0, 1, len(self._gRDF.index.levels[0].values)))
        self._plotWorker(self._gRDF                        ,\
                         self._scanNameDict                ,\
                         self._scanNormalizationOrUnitsDict,\
                         self._scanXLabel                  ,\
                         scanColors                        ,\
                        )

        # Make the mode plot
        gRDFSwap = self._gRDF.swaplevel()
        modeColors =\
            seqCMap2(np.linspace(0, 1, len(gRDFSwap.index.levels[0].values)))
        self._plotWorker(gRDFSwap                          ,\
                         self._modeNameDict                ,\
                         self._modeNormalizationOrUnitsDict,\
                         self._modeXLabel                  ,\
                         modeColors                        ,\
                        )
    #}}}

    #{{{_plotWorker
    def _plotWorker(self                    ,\
                    gRDF                    ,\
                    scanOrModeNameDict      ,\
                    normalizationOrUnitsDict,\
                    xlabel                  ,\
                    colors                  ,\
                    ):
        #{{{docstring
        """
        The worker which works for plotSaveGrowthRates.

        Parameters
        ----------
        gRDF : DataFrame
            The growth rates data frame.
        scanOrModeNameDict : dict
            Dictionary containing the plot name
        normalizationOrUnitsDict : dict
            Dictionary containing the normalization or the units
        xlabel : str
            String to use on the x-axis label.
        colors : array 2d
            Array describing the colors.
        """
        #}}}

        # Create the axes
        figSize = SizeMaker.array(1, 2, w = 6.3, aSingle = 0.7*0.7)
        fig, (imAx, reAx) =\
            plt.subplots(nrows=2, sharex=True, figsize=figSize)

        for outerInd, color in zip(gRDF.index.levels[0].values, colors):
            # Get the value for the legend
            pltOuterIndDict =\
                    {"val":plotNumberFormatter(outerInd, None)}
            label = self._legendTemplate.format(scanOrModeNameDict      ,\
                                                normalizationOrUnitsDict,\
                                                pltOuterIndDict)

            # Update the colors legend
            self._errorbarOptions.update({"color":color, "ecolor":color})

            yerr = None if self._analytic\
                    else gRDF.loc[outerInd]["growthRateStd"].values

            (_, caps, _) = imAx.errorbar(\
                gRDF.loc[outerInd]["growthRate"].index.values,\
                gRDF.loc[outerInd]["growthRate"].values      ,\
                yerr  = yerr                                 ,\
                label = label                                ,\
                **self._errorbarOptions)

            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            yerr = None if self._analytic\
                    else gRDF.loc[outerInd]["angularVelocityStd"].values

            (_, caps, _) = reAx.errorbar(\
                gRDF.loc[outerInd]["angularVelocity"].index.values,\
                gRDF.loc[outerInd]["angularVelocity"].values      ,\
                yerr  = yerr                                      ,\
                **self._errorbarOptions)

            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

        # Set labels
        imAx.set_ylabel(self._imLabel)
        reAx.set_ylabel(self._reLabel)

        # Set scan label
        reAx.set_xlabel(xlabel)

        # Set the title
        fig.suptitle(self._title)

        # Tweak the ticks
        # Add 10% margins for readability
        imAx.margins(x=0.1, y=0.1)
        reAx.margins(x=0.1, y=0.1)

        # Make the plot look nice
        self._ph.makePlotPretty(imAx, loc ="upper left")
        self._ph.makePlotPretty(reAx, loc ="upper left",\
                                legend = False, rotation = 45)

        # Set the ticks
        # reAx.tick_params(labelbottom="off")
        ticks = tuple(gRDF.loc[outerInd]["growthRate"].index.values)
        reAx.xaxis.set_ticks(ticks)
        imAx.xaxis.set_ticks(ticks)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.1, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            # Sets the save name
            fileName = "{}-{}.{}".\
                format(self._fileName, gRDF.index.names[0], self._extension)

            self._ph.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
