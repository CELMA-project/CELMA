#!/usr/bin/env python

"""
Gives blobs per unit time for the B-Scan.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper

sigMul = 3
scans = ("B0_0.06", "B0_0.08", "B0_0.1")
blobsPerTime = []
holesPerTime = []
scanVal = []
for scan in scans:
    path = "../CSDXMagFieldScanAr/visualizationPhysical/{}/blobs/".format(scan)
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
    curScan = float(scan[4:])
    blobs   = counts.blobs
    holes   = counts.holes
    blobsPerTime.append(blobs/lastT)
    holesPerTime.append(holes/lastT)
    scanVal.append(curScan)

    print("Found {} blobs and {} holes for {}".format(blobs, holes, scan))

fig, ax = plt.subplots(figsize = SizeMaker.standard(s=0.4))
ax.plot(scanVal, blobsPerTime, "k", ls="", alpha=0.7,\
        marker="o", ms=10, label="$\mathrm{Blobs}$")
ax.plot(scanVal, holesPerTime, "gray", ls="", alpha=0.7,\
        marker="d", ms=10, label="$\mathrm{Holes}$")
ax.set_ylabel(r"$\mathrm{Average\;counts\;per\;second}$")
ax.set_xlabel(r"$B [\mathrm{T}]$")
PlotHelper.makePlotPretty(ax, rotation = 45, loc = "upper right")

PlotHelper.savePlot(fig, "BScanBlobCount.pdf")
