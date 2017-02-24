#!/usr/bin/env python

"""
Tweak mode nr plot by changing x axis.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper
from CELMAPy.growthRates import PlotGrowthRates

path = "../../CSDXNeutralScanAr/visualizationPhysical/all/growthRates/n-growthRates-rho-0.0363-z-0.7-modeNr.pickle"
with open(path, "rb") as f:
    fig = pickle.load(f)

axes       = fig.get_axes()
imAx, reAx = axes

# Obtain the errobarOptions
errorbarOptions = PlotGrowthRates.errorbarOptions
errorbarOptions.pop("color" )
errorbarOptions.pop("ecolor")

# Used for transformation
n0 = 1e19

yLabels = []
for ax in axes:
    lines = ax.get_lines()
    meanLines  = lines[0::3]
    yDownLines = lines[1::3]
    yUpLines   = lines[2::3]

    x      = []
    meanY  = []
    yDownY = []
    yUpY   = []
    colors = []

    # Transform the x-axis, and cast is to per cents
    x = meanLines[0].get_xdata()
    x = n0/(n0+np.array(x))*100

    for mean, yDown, yUp in zip(meanLines, yDownLines, yUpLines):
        meanY .append(mean .get_ydata())
        yDownY.append(yDown.get_ydata())
        yUpY  .append(yUp  .get_ydata())

        colors.append(mean .get_color())

    # Clear axis, but save the label first
    yLabels.append(ax.get_ylabel())
    ax.cla()

    for y, yDown, yUp, color in zip(meanY, yDownY, yUpY, colors):
        yDown = np.array(yDown) - np.array(y)
        yUp   = np.array(y) - np.array(yUp)
        yerr  = (np.array(yDown), np.array(yUp))

        # Replot
        ax.errorbar(x                ,\
                    y                ,\
                    color = color    ,\
                    yerr  = yerr     ,\
                    **errorbarOptions)

        PlotHelper.makePlotPretty(ax, legend = False, rotation = 45)

# Remove x axis on imAx
imAx.set_xticklabels(imAx.get_xticklabels(), visible=False)

# Set proper ticks
reAx.xaxis.set_ticks(x)
reAx.set_xticklabels([r"${:d}\;\%$".format(int(i)) for i in x])

# Set the labels
imAx.set_ylabel(yLabels[0])
reAx.set_ylabel(yLabels[1])
reAx.set_xlabel(r"$\mathrm{Ionization\;degree}$")

PlotHelper.savePlot(fig, "growthRatesNnModes.pdf")
