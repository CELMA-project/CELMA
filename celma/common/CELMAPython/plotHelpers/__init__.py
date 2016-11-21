#!/usr/bin/env python

""" Init for the plotHelpers package """

from .plotHelper import PlotHelper
from .plotNumberFormatter import plotNumberFormatter
from .improvedCollect import collectiveCollect, safeCollect
from .derivatives import DDZ, DDX, findLargestRadialGrad
import matplotlib.pyplot as plt

# Set the plot style for all plots
titleSize = 30
# WARNING: This makes it slow, but is needed in order not to get
#          UserWarning: findfont: Font family ['serif'] not found
plt.rc("text",   usetex=True)

font = {"family":"serif","size":30, "serif": ["computer modern roman"]}
plt.rc("font",   **font)
plt.rc("axes",   labelsize = 25, titlesize = titleSize)
plt.rc("xtick",  labelsize = 25)
plt.rc("ytick",  labelsize = 25)
plt.rc("legend", fontsize  = 20)
plt.rc("lines",  linewidth = 2)

# Set proper backend
try:
    plt.figure(0)
except RuntimeError:
    plt.switch_backend("Agg")
    plt.figure(0)
plt.close(0)

# Set the colorfunc
seqCMap  = plt.get_cmap("inferno")
seqCMap2 = plt.get_cmap("plasma")
seqCMap3 = plt.get_cmap("viridis")
divCMap  = plt.get_cmap("BrBG")
