#!/usr/bin/env python

"""
Tweaks the statistics obtained from the blob runner.
"""

import pickle
import matplotlib.pylab as plt
import matplotlib.patches as pa
from matplotlib.ticker import MaxNLocator
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper, seqCMap2

# Set the label colors
colors = seqCMap2(np.linspace(0.25,0.75,3))

sigmas = (2,3,4)
sD = {key:{"waiting":None, "pulse":None, "color":c}\
            for key, c in zip(sigmas,colors)}

scan = "B0_0.08"
path = "../CSDXMagFieldScanAr/visualizationPhysical/{}/blobs/".format(scan)

# Obtain the patches
for key in sD.keys():
    fileName =\
        os.path.join(path, str(key), "temporalStats-blobs.pickle")
    with open(fileName, "rb") as f:
        fig = pickle.load(f)

    axes = fig.get_axes()
    wAx = axes[0]
    pAx = axes[1]
    wTitle = wAx.get_title()
    pTitle = pAx.get_title()

    wPatches = wAx.patches
    pPatches = pAx.patches

    title   = fig.texts[0].get_text()
    xLabel  = wAx.get_xlabel()
    yLabelP = pAx.get_ylabel()
    yLabelW = wAx.get_ylabel()

    sD[key]["waiting"] = wPatches
    sD[key]["pulse"]   = pPatches

# Make a new figure
fig, (wAx, pAx) = plt.subplots(ncols=2, figsize = SizeMaker.array(2,1))

patchesForLegend = []
for key in sD.keys():
    for nr, p in enumerate(sD[key]["waiting"]):
        # Make a new patch
        label = r"${}\sigma$".format(key) if nr == 0 else None
        wAx.add_patch(pa.Rectangle(p.get_xy(), p.get_width(),
                      p.get_height(), color = sD[key]["color"], alpha = 0.5,\
                      label = label))
    for p in sD[key]["pulse"]:
        # Make a new patch
        pAx.add_patch(pa.Rectangle(p.get_xy(), p.get_width(),
                      p.get_height(), color = sD[key]["color"], alpha =
                      0.5))

# Set the decorations
wAx.set_title(wTitle)
wAx.set_xlabel(xLabel)
wAx.set_ylabel(yLabelW)
pAx.set_title(pTitle)
pAx.set_xlabel(xLabel)
pAx.set_ylabel(yLabelP)

wAx.autoscale()
pAx.autoscale()
fig.suptitle(title, y=1.1)

PlotHelper.makePlotPretty(wAx, rotation = 45)
PlotHelper.makePlotPretty(pAx, rotation = 45, legend = None)

wAx.yaxis.set_major_locator(MaxNLocator(integer=True))
pAx.yaxis.set_major_locator(MaxNLocator(integer=True))

PlotHelper.savePlot(fig, "blobStats.pdf")
