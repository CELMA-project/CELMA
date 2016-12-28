#!/usr/bin/env python

"""
Contains class to plot the growth rates
"""

from ..plotHelpers import PlotHelper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import os











YOU ARE HERE
import pandas as pd
B0s=[0.01,0.02]
modes=[100,200,300]
multiTuples = []
fullDict = {"gr":[], "grStd":[], "aV":[], "aVStd":[]}
for B0 in B0s:
    for mode in modes:
        for key in fullDict.keys():
            fullDict[key].append(mode*B0)
        multiTuples.append((B0, mode))
df = pd.DataFrame(fullDict, index=pd.MultiIndex.from_tuples(multiTuples, names=["B0","modes"]))
print(df)
print("Outer level name")
print(df.index.names[0])
print("Get outer level indices")
print(df.index.levels[0])
# Loop through indices
# One of the indices is 0.01
import matplotlib.pyplot as plt
plt.errorbar(df.loc[0.01]["gr"].index, df.loc[0.01]["gr"].values, yerr=df.loc[0.01]["grStd"].values)
ax = plt.gca()
ax.set_ylabel("gr")
ax.set_xlabel(df.loc[0.01]["gr"].index.name)
#plt.show()
# Swap levels
df2 = df.swaplevel()
print("Get outer level name")
print(df2.index.names[0])




#{{{PlotGrowthRates
class PlotGrowthRates(PlotSuperClass):
    """Class which contains the growth rates and the plotting configuration."""

    #{{{Static members
    _mapToPltText = {"modeNr":"$\mathrm{Mode\quad number}$",\
                     "B0"    :"$B_0$"                      ,\
                     "Te0"   :"$T_e$"                      ,\
                     "nn"    :"$n_n$"                      ,\
                     "length":"$z$"                        ,\
                    }

    _errorbarOptions = {"color"     :"k",\
                        "fmt"       :"o",\
                        "fmt"       :"o",\
                        "markersize":10 ,\
                        "ecolor"    :"k",\
                        "capsize"   :7  ,\
                        "elinewidth":3  ,\
                        }

    _markeredgewidth = 3

    _all = slice(None)
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
        self._varName = self._varName.replace("Magnitude","")
        self._varName = self._varName.replace("AngularVelocity","")

        # Obtain the color
        self._colors = seqCMap2(np.linspace(0, 1, self._nModes))

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
        self._varLabelTemplate = r"${{}}${}".format(unitsOrNormalization)

        # Set the time label
        self._timeLabel = self._ph.tTxtDict["tTxtLabel"]
    #}}}

    #{{{plotSaveGrowthRates
    def plotSaveGrowthRates(self):
        """
        Performs the actual plotting.
        """

        keys = sorted(self._fourierModes.keys())

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
            for modeNr, color in zip(range(1, self._nModes), self._colors):
                label=r"$m_\theta={}$".format(modeNr)

                ax.plot(\
                    self._fourierModes[key]["time"],\
                    self._fourierModes[key][self._varName+"Magnitude"][:,modeNr],\
                    color=color, label=label)

            # Set axis labels
            ax.set_xlabel(self._timeLabel)
            ax.set_ylabel(self._varLabel)

            # Set logarithmic scale
            ax.set_yscale("log")

            # Set the title
            ax.set_title(title)

            # Make the plot look nice
            self._ph.makePlotPretty(ax, rotation = 45)

            if self._showPlot:
                plt.show()

            if self._savePlot:
                # Sets the save name
                fileName = "{}-rho-{}-z-{}.{}".\
                    format(self._fileName, rhoVal, zVal, self._extension)
                self._ph.savePlot(fig, fileName)

            plt.close(fig)
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
