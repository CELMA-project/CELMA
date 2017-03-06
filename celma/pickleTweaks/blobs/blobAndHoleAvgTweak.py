#!/usr/bin/env python

"""
Tweaks the averaged structures obtained from the blob runner.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper, seqCMap2

scan = "B0_0.08"
path = "../../CSDXMagFieldScanAr/visualizationPhysical/{}/blobs/".format(scan)

# Set the label colors
colors = seqCMap2(np.linspace(0.25,0.75,3))

sigmas = (2,3,4)
data = {key:{"bTime":None, "bDens":None, "hTime":None, "hDens":None}\
        for key in sigmas}

# Obtain the data
for key in data.keys():
    fileName =\
        os.path.join(path, str(key), "timeTrace-blobAndHole-avg-0.pickle")
    with open(fileName, "rb") as f:
        fig = pickle.load(f)

    axes = fig.get_axes()

    for nr, ax in enumerate(axes):
        blobOrHole = "b" if nr == 0 else "h"
        l = ax.get_lines()[0]
        time, dens = l.get_data()
        ax.lines.pop()
        data[key][blobOrHole + "Time"] = time
        data[key][blobOrHole + "Dens"] = dens

# Replot
axes = fig.get_axes()

for nr, ax in enumerate(axes):
    blobOrHole = "b" if nr == 0 else "h"
    for key, c in zip(data.keys(), colors):
        ax.plot(data[key][blobOrHole + "Time"],\
                data[key][blobOrHole + "Dens"],\
                color=c,
                label="${}\sigma$".format(key),\
                alpha=0.75)
        if nr == 0:
            leg = ax.legend(loc       = "best",\
                            fancybox  = True  ,\
                            numpoints = 1     ,\
                            )
            leg.get_frame().set_alpha(0.5)

PlotHelper.savePlot(fig, "blobsAndHoles-{}.pdf".format(scan))
