#!/usr/bin/env python

"""
Tweak scan plot by changing the legend.

NOTE: Couldn't remove the legend, so needed to make a new figure
"""

import pickle
import matplotlib.pylab as plt
import numpy as np
from matplotlib.pylab import MaxNLocator

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper
from CELMAPy.growthRates import PlotGrowthRates

path = "../../CSDXNeutralScanAr/visualizationPhysical/all/growthRates/n-growthRates-rho-0.0363-z-0.7-nn.pickle"
with open(path, "rb") as f:
    fig = pickle.load(f)

axes = fig.get_axes()
imAx, reAx = axes

# Used for transformation
n0 = 1e19

# Update the handles and labels, and recreate them
handles, labels = imAx.get_legend_handles_labels()

newLabels = []
for label in labels:
    nr = float(label.split("$")[3].replace("\cdot 10", "e").replace("^{","").replace("}",""))
    nr = n0/(n0+np.array(nr))*100
    newLabels.append(r"$d\;=\;{:d}\%$".format(int(np.ceil(nr))))

# Create new fig
newFig, newAxes = plt.subplots(nrows=2                      ,\
                               figsize=fig.get_size_inches(),\
                               sharex=True)

# Obtain the errobarOptions
errorbarOptions = PlotGrowthRates.errorbarOptions
errorbarOptions.pop("color" )
errorbarOptions.pop("ecolor")

for ax, newAx in zip(axes, newAxes):
    lines = ax.get_lines()
    meanLines  = lines[0::3]
    yDownLines = lines[1::3]
    yUpLines   = lines[2::3]

    x = meanLines[0].get_xdata()

    meanY  = []
    yDownY = []
    yUpY   = []
    colors = []

    for mean, yDown, yUp in zip(meanLines, yDownLines, yUpLines):
        meanY .append(mean .get_ydata())
        yDownY.append(yDown.get_ydata())
        yUpY  .append(yUp  .get_ydata())

        colors.append(mean .get_color())

    for y, yDown, yUp, color in zip(meanY, yDownY, yUpY, colors):
        yDown = np.array(yDown) - np.array(y)
        yUp   = np.array(y) - np.array(yUp)
        yerr  = (np.array(yDown), np.array(yUp))

        # Replot
        newAx.errorbar(x                ,\
                       y                ,\
                       color = color    ,\
                       yerr  = yerr     ,\
                       **errorbarOptions)

        newAx.set_ylabel(ax.get_ylabel())

newAxes[0].set_title (fig.texts[0].get_text())
newAxes[1].set_xlabel(reAx.get_xlabel())
PlotHelper.makePlotPretty(newAxes[0], legend = False, rotation = 45)
PlotHelper.makePlotPretty(newAxes[1], legend = False, rotation = 45)

# Manually set the x and the y as the figure has different coordinates
x = 0.4
y = -0.4

newFig.legend(handles,\
              newLabels,\
              bbox_to_anchor=(x, y),\
              ncol=2,\
              loc="upper center",\
              borderaxespad=0.,\
              bbox_transform = newAxes[1].transAxes,\
             )

newAxes[0].get_xaxis().set_major_locator(MaxNLocator(integer=True))

PlotHelper.savePlot(newFig, "growthRatesNnScan.pdf")
