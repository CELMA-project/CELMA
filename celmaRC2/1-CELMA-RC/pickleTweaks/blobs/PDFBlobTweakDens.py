#!/usr/bin/env python

"""
Tweaks the density PDF obtained from the blob runner.
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
path = "../../CSDXMagFieldScanAr/visualizationPhysical/{}/PDFs/".format(scan)
fileName = os.path.join(path, "PDF.pickle")
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

# Obtain the PDF stats
fileName = os.path.join(path, "PDFStats.pickle")
with open(fileName, "rb") as f:
    PDFStats = pickle.load(f)

posKey = list(PDFStats.keys())[0]
skew       = PDFStats[posKey]["skew"]
kurtExcess = PDFStats[posKey]["kurtExcess"]

# Add skewness and kurtosis to the plot
# Get the textSize
textPos = (0.725, 0.945)
# Add text
SKTxt = ax.\
        text(*textPos,\
            "$S = {:.2f}$\n $K = {:.2f}$".\
            format(skew, kurtExcess),\
            transform = ax.transAxes,\
            ha="left", va="top",\
            bbox={"facecolor":"white", "alpha":0.5, "pad":5},\
            )
SKTxt.set_fontsize(txtSize)

PlotHelper.savePlot(fig, "blobDensPDF_{}.pdf".format(scan))
