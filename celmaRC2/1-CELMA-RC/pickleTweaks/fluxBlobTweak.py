#!/usr/bin/env python

"""
Tweaks the flux obtained from the blob runner.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np
from scipy import stats

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper, seqCMap2

# Will not use namedTuple here as it is unnecessary complicated to
# pickle/unpickle
sigMul = (2,3,4)
PDFStats={"stdDev":{key:None for key in sigMul}, "skew":None, "kurtExcess":None}

scan = "B0_0.08"
path = "../CSDXMagFieldScanAr/visualizationPhysical/{}/blobs/".format(scan)
fileName = os.path.join(path, "3", "n-radialFluxes-fluct.pickle")
with open(fileName, "rb") as f:
    fig = pickle.load(f)

ax = fig.get_axes()[0]

# Make title from legend
handle, label = ax.get_legend_handles_labels()

fig.suptitle(label[0])

# Remove old legend
leg = ax.legend()
leg.remove()

# Recolor line
l = ax.get_lines()[0]
l.set_color("k")
time, flux = l.get_data()
t = (time[0], time[-1])

# Set the label colors
colors = seqCMap2(np.linspace(0.25,0.75,3))

# Calculate the standard deviation
stdDev = flux.std()

for mul, color in zip(sigMul, colors):
    fileName = os.path.join(path, str(mul), "counts.pickle")
    with open(fileName, "rb") as f:
        counts = pickle.load(f)
    curSig = stdDev*mul
    cond = (curSig,)*2
    label = r"${}\sigma$".format(mul)
    label += "\n" + "$\mathrm{{Blobs}}:{}$".format(counts.blobs)
    label += "\n" + "$\mathrm{{Holes}}:{}$".format(counts.holes)
    ax.plot(t,cond, alpha = 0.7, lw=2, label=label, color=color)
    PDFStats["stdDev"][mul] = curSig

# Pickle the PDFStats
PDFStats["skew"]       = stats.skew(flux)
PDFStats["kurtExcess"] = stats.kurtosis(flux)
fileName = os.path.join(path, "3", "PDFStats.pickle")
with open(fileName, "wb") as f:
    pickle.dump(PDFStats, f, pickle.HIGHEST_PROTOCOL)
    print("Pickled {}".format(fileName))

# Move legend outside
handles, labels = ax.get_legend_handles_labels()
# Remove handle and label deleted in the legend at the same time
handles = handles[1:]
labels  = labels [1:]
# Have them in reverse order for better plot
fig.legend(handles[::-1],\
           labels [::-1],\
           bbox_to_anchor=(1.05, 1.0),\
           loc="upper left",\
           borderaxespad=0.,\
           bbox_transform = ax.transAxes,\
           )

# Tweak the size a bit
fig.set_figwidth (4.0)
fig.set_figheight(2.5)

PlotHelper.savePlot(fig, "blobFluxTimeTrace_{}.pdf".format(scan))
