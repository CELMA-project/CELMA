#!/usr/bin/env python

"""
Tweak the skewness and kurtosis as a function of nn.
To be used in thesis.
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

from CELMAPy.plotHelpers import SizeMaker, PlotHelper, plotNumberFormatter, seqCMap3

scans  = (\
          "nn_9.9e+20"              ,\
          "nn_4e+19"                ,\
          "nn_1.5e+19"              ,\
          "nn_2.5e+18"              ,\
          "nn_6.666666666666668e+18",\
         )

markers = ("*", "o", "s", "v", "^")
ls = (\
      (0, ()),\
      (0, (5, 5)),\
      (0, (3, 1, 1, 1)),\
      (0, (3, 1, 1, 1, 1, 1)),\
      (0, (1, 1)),\
     )

colors = seqCMap3(np.linspace(0, 1, len(scans))[::-1])
sD = {s:{"ls":l,"marker":m,"color":c}\
        for s, m, l, c in zip(scans, markers, ls, colors)}

for scan in scans:
    folder = "../../CSDXNeutralScanAr/visualizationPhysical/{}/skewKurt".format(scan)
    picklePath = os.path.join(folder, "skewKurt.pickle")
    with open(picklePath, "rb") as f:
        fig = pickle.load(f)

    ax = fig.get_axes()[0]

    # Get the legends, and make keys out of these
    handles, labels = ax.get_legend_handles_labels()

    # Obtain the lines
    for nr, label in enumerate(labels):
        sD[scan][label] = ax.get_lines()[nr].get_data()[1]

    xAxis = ax.get_lines()[0].get_data()[0]

    plt.close(fig)

# Make a new figure
fig, (sAx, kAx) = plt.subplots(nrows=2, figsize = SizeMaker.standard(w=3, a=1.5),\
                               sharex=True)

for scan in scans:
    curScan = float(scan[3:])
    label = r"$n_n =$ {}".format(plotNumberFormatter(curScan, None)) +\
            r" $\mathrm{m}^{-3}$"
    for key in sD[scan].keys():
        if "skew" in key.lower():
            sAx.plot(xAxis, sD[scan][key],\
                     ls     = sD[scan]["ls"],\
                     color  = sD[scan]["color"],\
                     marker = sD[scan]["marker"],\
                     ms     = 7,\
                     alpha  = 0.7,\
                     label  = label\
                     )
            sAx.set_ylabel(key.replace("Skewness", r"Skewness \;"))
        elif "kurt" in key.lower():
            kAx.plot(xAxis, sD[scan][key],\
                     ls     = sD[scan]["ls"],\
                     color  = sD[scan]["color"],\
                     marker = sD[scan]["marker"],\
                     alpha  = 0.7,\
                     ms     = 7,\
                     )
            kAx.set_ylabel(key.replace("quadkurtosis", r"; kurtosis \;"))
            kAx.set_xlabel(r"$\rho$ $[m]$")

PlotHelper.makePlotPretty(sAx, rotation = 45)
PlotHelper.makePlotPretty(kAx, rotation = 45, legend=None)

sAx.legend(bbox_to_anchor=(1.6, 0.5),\
           loc="center",\
           borderaxespad=0.,\
           bbox_transform = sAx.transAxes,\
           )

fig.tight_layout()

fileName = "nnScanSkewKurt.pdf"
PlotHelper.savePlot(fig, fileName)
