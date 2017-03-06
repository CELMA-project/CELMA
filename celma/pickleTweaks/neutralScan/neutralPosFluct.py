#!/usr/bin/env python

"""
Tweak of profile plots as a function of nn.
"""

import pickle
import matplotlib.pylab as plt
import matplotlib.lines as mlines
import numpy as np
from subprocess import Popen

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper, seqCMap3

scans  = (\
          "nn_4e+19"                ,\
          "nn_1.5e+19"              ,\
          "nn_6.666666666666668e+18",\
          "nn_2.5e+18"              ,\
         )
markers = ("*", "o", "s", "v", "^")
ls = (\
      (0, ()),\
      (0, (5, 5)),\
      (0, (3, 1, 1, 1)),\
      (0, (3, 1, 1, 1, 1, 1)),\
      (0, (1, 1)),\
     )

colors = seqCMap3(np.linspace(2/5, 1, len(scans))[::-1])
sD = {s:{"ls":l,"marker":m,"color":c}\
        for s, m, l, c in zip(scans, markers, ls, colors)}

for scan in scans:
    folder = "../../CSDXNeutralScanAr/visualizationPhysical/{}/radialProfiles".format(scan)
    picklePath = os.path.join(folder, "n-phi-posOfFluct.pickle")
    with open(picklePath, "rb") as f:
        fig = pickle.load(f)

    ax = fig.get_axes()[2]

    # Get the line
    sD[scan]["line"] = ax.get_lines()[0].get_data()[1]
    xAxis = ax.get_lines()[0].get_data()[0]

    # Set legend to ylabel
    handles, labels = ax.get_legend_handles_labels()
    yLabel = labels[0]
    xLabel = fig.get_axes()[3].get_xlabel()

    plt.close(fig)

# Make a new figure
fig, ax = plt.subplots(figsize = SizeMaker.standard(s=0.5))

# Used for conversion
n0 = 1e19

for scan in scans:
    curScan = n0/(n0+float(scan[3:]))*100
    ax.plot(xAxis, sD[scan]["line"],\
            ls     = sD[scan]["ls"],\
            color  = sD[scan]["color"],\
            marker = sD[scan]["marker"],\
            ms     = 7,\
            alpha  = 0.7,\
            label  = r"$d = {:d} \;\%$".\
                     format(int(curScan))\
            )

ax.set_ylabel(yLabel)
ax.set_xlabel(xLabel)

PlotHelper.makePlotPretty(ax, rotation = 45)

fileName = "nnScanPosOfFluct.pdf"
PlotHelper.savePlot(fig, fileName)
