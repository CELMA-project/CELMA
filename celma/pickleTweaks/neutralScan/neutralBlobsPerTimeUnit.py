#!/usr/bin/env python

"""
Gives blobs per unit time for the nn-Scan.
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

sigMul = 3
scans  = (\
          "nn_4e+19"                ,\
          "nn_1.5e+19"              ,\
          "nn_6.666666666666668e+18",\
          "nn_2.5e+18"              ,\
         )
blobsPerTime = []
holesPerTime = []
scanVal = []
for scan in scans:
    path = "../../CSDXNeutralScanAr/visualizationPhysical/{}/blobs/".format(scan)
    # Obtain the counts
    fileName = os.path.join(path, str(sigMul), "counts.pickle")
    with open(fileName, "rb") as f:
        counts = pickle.load(f)

    # Obtain the time
    fileName = os.path.join(path, str(sigMul), "n-radialFluxes-fluct.pickle")
    with open(fileName, "rb") as f:
        fig = pickle.load(f)

    ax = fig.get_axes()[0]
    time, _ = ax.get_lines()[0].get_data()
    lastT = time[-1]
    curScan = float(scan[3:])
    blobs   = counts.blobs
    holes   = counts.holes
    blobsPerTime.append(blobs/lastT)
    holesPerTime.append(holes/lastT)
    scanVal.append(curScan)

    print("Found {} blobs and {} holes for {}".format(blobs, holes, scan))

plt.close("all")

fig, ax = plt.subplots(figsize = (4,1))

# Format scan vals
n0      = 1e19
scanVal = n0/(n0+np.array(scanVal))*100

yTicks    = tuple(r"${:d} \; \%$".format(int(d)) for d in scanVal)
xBarVals  = tuple(range(len(blobsPerTime)*3))
blobsVals = xBarVals[0::3]
holesVals = xBarVals[1::3]
spaceVals = xBarVals[2::3]

# Find the lowest value
lowestVal = np.min(holesPerTime if np.min(holesPerTime) < np.min(blobsPerTime)\
                                else blobsPerTime)
# Pad with 5%
lowestVal = 0.95*lowestVal

# Get the color cycler
# http://stackoverflow.com/questions/13831549/get-matplotlib-color-cycle-state
prop_cycler = ax._get_lines.prop_cycler
# Plot
blobRects = ax.barh(blobsVals       ,\
                    blobsPerTime    ,\
                    left = lowestVal,\
                    color="black"   ,\
                   )

holeRects = ax.barh(holesVals        ,\
                    holesPerTime     ,\
                    left = lowestVal ,\
                    color="white"    ,\
                    edgecolor="black",\
                   )

# Add space
ax.barh(spaceVals       ,\
        blobsPerTime    ,\
        left = lowestVal,\
        color="white"   ,\
       )

# Set text
for nr, (rect, txt) in enumerate(zip(holeRects, yTicks)):
    width  = rect.get_width()
    height = rect.get_y()
    # Add 4% padding
    x      = 0
    y      = height
    ax.text(x, y, txt, ha="center", va="center")

ax.set_xlabel(r"$[\mathrm{s}^{-1}]$")

ax.spines["top"]  .set_visible(False)
ax.spines["left"] .set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_yaxis().set_visible(False)

PlotHelper.makePlotPretty(ax, legend=False, rotation=45)
ax.grid(False)

# Make legend manually
handles = (holeRects[0], blobRects[0])
labels  = ("$\mathrm{Holes}$",\
           "$\mathrm{Blobs}$",\
           )

fig.suptitle(r"$\mathrm{Average\;blobs\;and\;holes\;per\;second}$", y=1.1)
fig.legend(handles                      ,\
           labels                       ,\
           bbox_to_anchor=(1.02, 1.0)   ,\
           loc="upper left"             ,\
           borderaxespad=0.             ,\
           bbox_transform = ax.transAxes,\
           )

PlotHelper.savePlot(fig, "nnScanBlobCount.pdf")
