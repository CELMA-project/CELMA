#!/usr/bin/env python

"""
Makes the plots to be used in the matrix showing the averaged blob.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np
from glob import glob
from subprocess import Popen

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper

fluct = "-fluct" # Set to "" if no fluctuation
blob = "blobs"
mode = "perp"
scan = "B0_0.06"
condition = "3"
path = "../../CSDXMagFieldScanAr/visualizationPhysical/{}/blobs/{}/".format(scan, condition)

# Obtain the marker
print("    Obatining marker position")
with open(os.path.join(path, "n-radialFluxes-fluct.pickle"), "rb") as curF:
    fig = pickle.load(curF)
    ax = fig.get_axes()[0]
    _, legend = ax.get_legend_handles_labels()
    rho   = float(legend[0].replace("$","").split("=")[1].split(" ")[1])
    theta = float(legend[0].replace("$","").split("=")[2].split(" ")[0].\
            replace("^{\\circ},", ""))*np.pi/180
    x     = rho*np.cos(theta)
    y     = rho*np.sin(theta)

    plt.close(fig)

# Find files matching the criteria
searchFor =\
    os.path.join(path, "n-{}-2D{}-{}-avg-*.pickle".format(mode, fluct, blob))

files = sorted(glob(searchFor))
if len(files) != 9:
    raise RuntimeError("Need the number of files to be exactly 9")

# Make a folder for the results
endFolder = "matrix-{}-{}-{}-{}".format(mode, blob, scan, fluct[1:])
if not(os.path.exists(endFolder)):
    os.makedirs(endFolder)

# Loop through the files
for nr, f in enumerate(files):
    print("    Opening {}".format(f))
    with open(f, "rb") as curF:
        fig = pickle.load(curF)

    ax, cbar = fig.get_axes()
    # Remove the colorbar
    fig.delaxes(cbar)

    # Keep the y-ticks only on figure 0, 3, 6
    if not(nr%3 == 0):
        ax.set_yticklabels([])
        ax.set_ylabel("")

    # Keep the x-ticks only on figure 6, 7, 8
    if not(nr >= 6):
        ax.set_xticklabels([])
        ax.set_xlabel("")

    # Redo the axis title
    t = r"$\tau={}$".format(ax.get_title().split("=")[-1][:-1].replace("$",""))
    ax.set_title(t.replace(r"\cdot 10^{-3}", r"\mu"))
    fig.set_size_inches(2.1, 2.1)

    # Add the marker
    ax.plot(x,y,"ko", markersize = 5, alpha=0.5)

    fileName = os.path.join(endFolder, "{}.pdf".format(nr))
    PlotHelper.savePlot(fig, fileName)
    Popen("pdfcrop {0} {0}".format(fileName), shell=True).wait()

# Save the colorbar
with open(f, "rb") as curF:
    fig = pickle.load(curF)

axes = fig.get_axes()
# Remove the plot
fig.delaxes(axes[0])
fileName = os.path.join(endFolder, "colorbar.pdf")
PlotHelper.savePlot(fig, fileName)
Popen("pdfcrop {0} {0}".format(fileName), shell=True).wait()
