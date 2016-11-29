#!/usr/bin/env python

""" Init for the plotHelpers package """

from .plotHelper import PlotHelper
from .plotNumberFormatter import plotNumberFormatter
from .improvedCollect import collectiveCollect, safeCollect
from .derivatives import DDZ, DDX, findLargestRadialGrad
import matplotlib.pyplot as plt
import os

# Set the plot style for all plots
titleSize = 30
size = 30

plt.rc("font",   size = size)
plt.rc("axes",   labelsize = 25, titlesize = titleSize)
plt.rc("xtick",  labelsize = 25)
plt.rc("ytick",  labelsize = 25)
plt.rc("legend", fontsize  = 20)
plt.rc("lines",  linewidth = 2)

# Set proper backend
try:
    plt.figure()
except RuntimeError:
    plt.switch_backend("Agg")
    plt.figure()

oldFont = {"family":plt.rcParams["font.family"],\
           "serif":plt.rcParams["font.serif"]}
try:
    # WARNING: This makes it slow, but is needed in order not to get
    #          UserWarning: findfont: Font family ['serif'] not found
    # Requires LaTeX and dvipng and Ghostscript installed
    plt.rc("text", usetex=True)

    font = {"family":"serif", "serif": ["computer modern roman"]}
    plt.rc("font", **font)

    fig, ax = plt.subplots()
    fileName = "tmp.png"
    plt.savefig(fileName)
    os.remove(fileName)
except (RuntimeError, FileNotFoundError) as er:
    plt.rc("text", usetex=False)
    plt.rc("font", **oldFont)

# Close all plots
plt.close("all")

# Set the colorfunc
seqCMap  = plt.get_cmap("inferno")
seqCMap2 = plt.get_cmap("plasma")
seqCMap3 = plt.get_cmap("viridis")
divCMap  = plt.get_cmap("BrBG")
