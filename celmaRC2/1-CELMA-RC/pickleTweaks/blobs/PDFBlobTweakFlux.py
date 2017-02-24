#!/usr/bin/env python

"""
Tweaks the flux PDF obtained from the blob runner.

To be runned after fluxBlobTweak.py as needs the numbers from the sigmas
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
fileName = os.path.join(path, "3", "PDF.pickle")
with open(fileName, "rb") as f:
    fig = pickle.load(f)

ax = fig.get_axes()[0]

# Make title from legend
handle, label = ax.get_legend_handles_labels()

fig.suptitle(label[0])

# Remove old legend
leg = ax.legend()
txtSize = leg.get_texts()[0].get_fontsize()
leg.remove()

# Recolor line
l = ax.get_lines()[0]
l.set_color("k")
time, PDFLine = l.get_data()
PDFLim = ax.get_ylim()

# Set the label colors
colors = seqCMap2(np.linspace(0.25,0.75,3))

# Obtain the PDF Stats
fileName = os.path.join(path, "3", "PDFStats.pickle")
with open(fileName, "rb") as f:
    PDFStats = pickle.load(f)

# Add lines to axis
stdDevKeys = sorted(PDFStats["stdDev"].keys())
for stdKey, color in zip(stdDevKeys, colors):
    std = (PDFStats["stdDev"][stdKey],)*2
    label = r"${}\sigma$".format(stdKey)
    ax.plot(std, PDFLim, alpha = 0.7, lw=2, label=label, color=color)

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

# Fix labels
old = r"u_{E\times B,\rho}}"
new = r"}\widetilde{u}_{E\times B,\rho}"
ylabel = ax.get_ylabel()
ylabel = ylabel.replace(old, new)
ax.set_ylabel(ylabel)
xlabel = ax.get_xlabel()
xlabel = xlabel.replace(old, new)
ax.set_xlabel(xlabel)

# Add skewness and kurtosis text
# Get the textSize
textPos = (0.725, 0.945)
# Add text
SKTxt = ax.\
        text(*textPos,\
            "$S = {:.2f}$\n $K = {:.2f}$".\
            format(PDFStats["skew"], PDFStats["kurtExcess"]),\
            transform = ax.transAxes,\
            ha="left", va="top",\
            bbox={"facecolor":"white", "alpha":0.5, "pad":5},\
            )
SKTxt.set_fontsize(txtSize)

PlotHelper.savePlot(fig, "blobFluxPDF_{}.pdf".format(scan))
