#!/usr/bin/env python

""" Init for the plotHelpers package """

from .plotHelper import PlotHelper
from .plotNumberFormatter import plotNumberFormatter
from .maxMinHelper import (getMaxMinAnimation,\
                           getLevelsAnimation,\
                           getVmaxVminLevels)
from .sizeMaker import SizeMaker
import os
import matplotlib.pyplot as plt

# Set proper backend from display
try:
    os.environ["DISPLAY"]
except KeyError:
    plt.switch_backend("Agg")

# Set the plot style for all plots
titleSize = 14

plt.rc("axes",   labelsize  = 12, titlesize = titleSize)
plt.rc("xtick",  labelsize  = 12)
plt.rc("ytick",  labelsize  = 12)
plt.rc("legend", fontsize   = 12)
plt.rc("lines",  linewidth  = 1)
plt.rc("lines",  markersize = 1.5)

# Use computer modern as default font
plt.rc("mathtext", fontset = "cm")

# Set the colorfunc
seqCMap   = plt.get_cmap("inferno")
seqCMap2  = plt.get_cmap("plasma")
seqCMap3  = plt.get_cmap("viridis")
divCMap   = plt.get_cmap("BrBG")
qualCMap  = plt.get_cmap("Vega10")
