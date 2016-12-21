#!/usr/bin/env python

"""
Contains class to plot the growth rates
"""

from ..plotHelpers import PlotHelper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import os

#{{{PlotGrowthRates
class PlotGrowthRates(object):
    """Class which contains the growth rates and the plotting configuration."""

    #{{{__init___
    def __init__(self                      ,\
                 paths                     ,\
                 growthRates               ,\
                 convertToPhysical = False ,\
                 showPlot          = False ,\
                 savePlot          = False ,\
                 extension         = "png" ,\
                 savePath          = "."   ,\
                 pltSize           = (9,12),\
                 ):
        #{{{docstring
        """
        The constructor for the PlotGrowthRates object.

        Sets the member data.

        Parameters
        ----------
        paths : str
            Paths to collect from (used to make the PlotHelper object)
        growthRates : DataFrame
            The data frame to be plotted. The data frame must be
            orgnized such that:
            1. The x axis will be given by splitting growthRates.columns
               by "=" and take the value right og the "="
            2. The index only has two levels
            3. The plot name is given at the first level of
               growthRates.index
            4. The growth rate, angular frequency and their standar
               deviation is given in the innermost level of growthRates
        showPlot : bool
            If the plots should be displayed.
        savePlot : bool
            If the plots should be saved.
        extension : str
            Extension to use on the plots
        savePath : str
            Path to save destination. Must exist.
        pltSize : tuple
            Size of the plots given as (x, y)
        """
        #}}}

        # Set the member data
        self._df        = growthRates
        self._showPlot  = showPlot
        self._savePlot  = savePlot
        self._extension = extension
        self._savePath  = savePath

        # Make the PlotHelper object
        self._helper = PlotHelper(paths[0][0]                          ,\
                                  useSpatial        = False            ,\
                                  convertToPhysical = convertToPhysical,\
                                 )

        # Set the plot size
        self._pltSize = pltSize

        # Dict used as an equivalent to C++ maps
        self._mapToPltText = {"modeNr":"$\mathrm{Mode\quad number}$",\
                              "B0"    :"$B_0$"                      ,\
                              "Te0"   :"$T_e$"                      ,\
                              "nn"    :"$n_n$"                      ,\
                              "length":"$z$"                        ,\
                             }

        self._errorbarOptions = {"color"     :"k",\
                                 "fmt"       :"o",\
                                 "fmt"       :"o",\
                                 "markersize":10 ,\
                                 "ecolor"    :"k",\
                                 "capsize"   :7  ,\
                                 "elinewidth":3  ,\
                                 }

        self._markeredgewidth = 3
        self._all = slice(None)
    #}}}

    #{{{plotGrowthRates
    def plotGrowthRates(self):
        #{{{docstring
        """
        Plots the growth rates.
        """
        #}}}

        # Loop through the figures
        for plotLabel in self._df.index.levels[0]:
            fig = plt.figure(figsize = self._pltSize)
            gs  = GridSpec(nrows=2, ncols=1)

            imAx   = fig.add_subplot(gs[0])
            realAx = fig.add_subplot(gs[1], sharex=imAx)

            plotLabelSplit = plotLabel.split("=")

            # We can now access unsing the loc method
            # The syntax is
            # self._df.loc[("indexLevel1", "indexLevel2", ...),\
            #              ("columnsLevel1", "columnsLevel2", ...)]
            # http://pandas.pydata.org/pandas-docs/stable/cookbook.html#cookbook-selection
            try:
                xAxis = tuple(float(txt.split("=")[1])\
                        for txt in\
                        self._df.loc[(plotLabel, "growthRate"), self._all].\
                        index)
            except KeyError as ke:
                message="{0}{1}WARNING: Only NaNs found in {2}. Skipping{1}{0}"
                print(message.format("\n", "!"*4, ke.args[0]))
                continue

            yAxisIm=self._df.loc[(plotLabel,"growthRate"), self._all].values
            yErrIm =self._df.loc[(plotLabel,"growthRateStd"), self._all].values
            yAxisRe=self._df.loc[(plotLabel,"angFreq"), self._all].values
            yErrRe =self._df.loc[(plotLabel,"angFreqStd"), self._all].values

            indexTxt =\
                self._df.loc[(plotLabel, "growthRate"), self._all].index[0].\
                                                                split("=")[0]

            # Convert units
            plotLabelSplit[-1], plotLabelNorm, plotLabelUnits =\
                self._helper.physicalUnitsConverter(plotLabelSplit[-1],\
                                                    plotLabelSplit[0])
            # NOTE: The whole dataframe can be multiplied using
            #       df.multiply, but we here only multiply the array for
            #       cleaner "helper" code
            xAxis, xAxisNorm, xAxisUnits =\
                self._helper.physicalUnitsConverter(xAxis  , indexTxt)
            yAxisIm, yAxisNorm, yAxisUnits =\
                self._helper.physicalUnitsConverter(yAxisIm, "growthRate")
            yErrIm, _, _  =\
                self._helper.physicalUnitsConverter(yErrIm , "growthRate")
            yAxisRe, _, _ =\
                self._helper.physicalUnitsConverter(yAxisRe, "growthRate")
            yErrRe, _, _  =\
                self._helper.physicalUnitsConverter(yErrRe , "growthRate")

            # Plot the growth rates
            (_, caps, _) = imAx.errorbar(xAxis,\
                                         yAxisIm,\
                                         yerr=yErrIm,\
                                         **self._errorbarOptions)

            # Set capsize
            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            # Angular frequencies
            (_, caps, _) = realAx.errorbar(xAxis,\
                                           yAxisRe,\
                                           yerr=yErrRe,\
                                           **self._errorbarOptions)

            # Set capsize
            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            # Add 10% margins for readability
            imAx  .margins(x=0.1, y=0.1)
            realAx.margins(x=0.1, y=0.1)

            rot = None if indexTxt == "modeNr" else 45
            PlotHelper.makePlotPretty(imAx,   yprune = "both", rotation = rot)
            PlotHelper.makePlotPretty(realAx, yprune = "both", rotation = rot)

            # Set the text
            if self._helper.convertToPhysical:
                if plotLabelSplit[0] == "modeNr":
                    suptitle = "{}$={}$".\
                        format(self._mapToPltText[plotLabelSplit[0]],\
                               plotLabelSplit[1])
                else:
                    suptitle = "{}$={}$ $[{}]$".\
                        format(self._mapToPltText[plotLabelSplit[0]],\
                               plotLabelSplit[1],\
                               plotLabelUnits)
                imLabel = r"$\omega_I$ $[{}]$".format(yAxisUnits)
                reLabel = r"$\omega_R$ $[{}]$".format(yAxisUnits)
            else:
                suptitle = "{}$={}{}$".\
                    format(self._mapToPltText[plotLabelSplit[0]],\
                           plotLabelSplit[1],\
                           plotLabelNorm)
                imLabel = r"$\omega_I{}$".format(yAxisNorm)
                reLabel = r"$\omega_R{}$".format(yAxisNorm)

            fig.suptitle(suptitle)
            imAx  .set_ylabel(imLabel)
            realAx.set_ylabel(reLabel)

            imAx.tick_params(labelbottom="off")

            # Set xlabel
            if indexTxt == "modeNr":
                xlabel = self._mapToPltText[indexTxt]
            else:
                if self._helper.convertToPhysical:
                    xlabel = "{} $[{}]$".\
                            format(self._mapToPltText[indexTxt], xAxisUnits)
                else:
                    xlabel = "{}{}".\
                            format(self._mapToPltText[indexTxt], xAxisNorm)

            realAx.set_xlabel(xlabel)

            # Adjust the subplots
            fig.subplots_adjust(hspace=0)

            # Sort the xAxis and yAxis, remove NaN's, and use them as ticks
            xAxis, yAxisIm = zip(*sorted(zip(xAxis, yAxisIm)))
            nonNan = np.where(np.isfinite(yAxisIm))[0]
            # Cut the NaN values from the xAxis (plus 1 as slice
            # excludes the last)
            ticks = tuple(xAxis[nonNan[0]:nonNan[-1]+1])

            # Set the ticks
            realAx.xaxis.set_ticks(ticks)
            imAx  .xaxis.set_ticks(ticks)

            if self._savePlot:
                fileName = "{}.{}".\
                    format(os.path.join(self._savePath,\
                           "growthrate_{}_{}".format(plotLabel, indexTxt)),\
                           self._extension)
                self._helper.savePlot(fig, fileName)

            if self._showPlot:
                plt.show()
    #}}}
#}}}
