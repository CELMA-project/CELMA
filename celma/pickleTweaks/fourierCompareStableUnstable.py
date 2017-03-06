#!/usr/bin/env python

"""
Tweak of Fourier transformed time traces used in Boussinesq compare.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper

picklePath = "../CSDXMagFieldScanAr/visualizationPhysical/B0_0.1/fourierModes/n-fourierModes-rho-0.0363-z-0.7.pickle"
with open(picklePath, "rb") as f:
    fig = pickle.load(f)

_, height = fig.get_size_inches()
width   = 3
fig.set_size_inches(width, height)
fileName = "FFT01.pdf"
PlotHelper.savePlot(fig, fileName)
plt.close(fig)

picklePath = "../CSDXMagFieldScanAr/visualizationPhysical/B0_0.02/fourierModes/n-fourierModes-rho-0.0388-z-0.7.pickle"
with open(picklePath, "rb") as f:
    fig = pickle.load(f)

fig.set_size_inches(width, height)
fileName = "FFT002.pdf"
PlotHelper.savePlot(fig, fileName)
plt.close(fig)
