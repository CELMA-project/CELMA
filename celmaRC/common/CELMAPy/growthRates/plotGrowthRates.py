#!/usr/bin/env python

"""
Contains class to plot the growth rates
"""

from ..superClasses import PlotSuperClass, seqCMap2, seqCMap3
from ..plotHelpers import plotNumberFormatter
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
    def setData(self, growthRateDataFrame, positionTuple):
        #{{{docstring
        """
        Sets the data to be plotted.

        This function:
            * Sets the data frame
            * Sets the variable labels
            * Prepares the save name

        Parameters
        ----------
        growthRateDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "growthRate"
                * "growthRateStd"
                * "averageAngularVelocity"
                * "averageAngularVelocityStd"
            over the observation "modeNr" over the observation "Scan"
        positionTuple : tuple
            The tuple containing (rho, z).
        """
        #}}}

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
        self._title(r"{}$,$ {}".format(rhoTitle, zTitle))

        # Set the growth rate
        self._gRDF = growthRateDataFrame

        # The variable name is stored in the outermost level of the data
        # frame observations
        self._varName = growthRateDataFrame.index.names[0]

        # Set the plot name
        self._pltVarName = self._ph.getVarPltName(self._varName)

        self._prepareLabels()

        # Set the var label
        self._varLabel = self._varLabelTemplate.\
            format(self._pltVarName, **self.uc.conversionDict[self._varName])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-rho-{}-z-{}".\
                format(self._varName, "growthRates", rhoVal, zVal))

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
            legendTemplate = r"${{scanOrMode}} = {{{{val}}}}{units}$"
        else:
            unitsOrNormalization = "${normalization}$"
            legendTemplate = r"${{scanOrMode}}{normalization} = {{{{val}}}}$"
        self._varLabelTemplate = r"${{}}${}".format(unitsOrNormalization)

        # Set the y-axes
        imag = r"$\Im(\omega_{{}})$".format(self._pltVarName.replace("$",""))
        real = r"$\Re(\omega_{{}})$".format(self._pltVarName.replace("$",""))
        self._imLabel =\
            self._varLabelTemplate.\
                format(imag, **self.uc.conversionDict["growthRate"])
        self._reLabel =\
            self._varLabelTemplate.\
                format(real, **self.uc.conversionDict["growthRate"])

        # Get the plot decoration for the scan
        scan            = self._gRDF.index.names[0]
        scanPltName     = self._ph.getVarPltName(scan)
        scanPltNameDict = {"scanOrMode":scanPltName}
        self._scanLegendTemplate = legendTemplate.\
            format(**self.uc.conversionDict[scan]).format(**scanPltNameDict)

        # Get the plot decoration for the modes
        mode            = self._gRDF.index.names[1]
        modePltName     = self._ph.getVarPltName(mode)
        modePltNameDict = {"scanOrMode":modePltName}
        self._modeLegendTemplate = legendTemplate.\
            format(**self.uc.conversionDict[mode]).format(**modePltNameDict)

        # Get the x labels
        self._scanXLabel = self._varLabelTemplate.\
                format(modePltName, **self.uc.conversionDict[mode])
        self._modeXLabel = self._varLabelTemplate.\
                format(scanPltName, **self.uc.conversionDict[scan])
    #}}}

    #{{{plotSaveGrowthRates
    def plotSaveGrowthRates(self):
        """
        Performs the plotting by calling the workers.
        """

        # Make the scan plot
        scanColors = seqCMap3(np.linspace(0, 1, len(self._gRDF.index.levels)))
        self._plotWorker(self._gRDF              ,\
                         self._scanLegendTemplate,\
                         self._scanXLabel        ,\
                         scanColors              ,\
                        )

        # Make the mode plot
        gRDFSwap = self._gRDF.swaplevel()
        modeColors = seqCMap2(np.linspace(0, 1, len(gRDFSwap.index.levels)))
        self._plotWorker(gRDFSwap                ,\
                         self._modeLegendTemplate,\
                         self._modeXLabel        ,\
                         modeColors              ,\
                        )
    #}}}

    #{{{_plotWorker
    def _plotWorker(self          ,\
                    gRDF          ,\
                    legendTemplate,\
                    xlabel        ,\
                    colors        ,\
                    ):
        #{{{docstring
        """
        The worker which works for plotSaveGrowthRates.

        Parameters
        ----------
        gRDF : DataFrame
            The growth rates data frame.
        legendTemplate : str
            String missing the keyword "val".
        xlabel : str
            String to use on the x-axis label.
        colors : array 2d
            Array describing the colors.
        """
        #}}}

        # Create the axes
        fig, (imAx, reAx) = plt.subplots(nrows=2, sharex=True)

        for outerInd, color in zip(gRDF.index.levels, colors):
            # Get the value for the legend
            pltOuterIndDict = {"val":plotNumberFormatter(outerInd, None)}

            # Update the colors legend
            self._errorbarOptions.update({"color":color, "ecolor":color})

            (_, caps, _) = imAx.errorbar(\
                    gRDF.loc[outerInd]["growthRate"].index            ,\
                    gRDF.loc[outerInd]["growthRate"].values           ,\
                    yerr  = gRDF.loc[outerInd]["growthRateStd"].values,\
                    label = legendTemplate.format(**pltOuterIndDict)  ,\
                    **self._errorbarOptions)

            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            (_, caps, _) = reAx.errorbar(\
                    gRDF.loc[outerInd]["averageAngularVelocity"].index       ,\
                    gRDF.loc[outerInd]["averageAngularVelocity"].values      ,\
                    yerr  =\
                       gRDF.loc[outerInd]["averageAngularVelocityStd"].values,\
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

        reAx.tick_params(labelbottom="off")
        ticks = tuple(gRDF.loc[outerInd]["growthRate"].index)

        # Set the ticks
        reAx.xaxis.set_ticks(ticks)
        imAx.xaxis.set_ticks(ticks)

        # Make the plot look nice
        self._ph.makePlotPretty(imAx)
        self._ph.makePlotPretty(reAx, rotation = 45)

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
