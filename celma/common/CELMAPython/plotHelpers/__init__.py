#!/usr/bin/env python

""" Init for the plotHelpers package """

from .plotHelper import PlotHelper
from .plotNumberFormatter import plotNumberFormatter
from .collectiveCollect import collectiveCollect
from .derivatives import DDZ, DDX, findLargestRadialGrad
import matplotlib.pyplot as plt

# Set the plot style for all plots
titleSize = 30
plt.rc("font",   size      = 30)
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
seqCMap  = plt.get_cmap("viridis")
seqCMap2 = plt.get_cmap("plasma")
divCMap  = plt.get_cmap("BrBG")
