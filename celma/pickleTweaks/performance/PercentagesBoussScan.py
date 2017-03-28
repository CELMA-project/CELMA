#!/usr/bin/env python

"""
Plots the percentages for the Boussinesq-scan.
"""

import pickle
import matplotlib.pylab as plt

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper

B0     = 0.08
modes  = ("Init", "Expand", "Linear", "Turbulence")

# Make place holders for the data
keys = ("Arithmetics"       ,\
        "Communication"     ,\
        "Input/output"      ,\
        "Laplace inversions",\
        "Time solver"       )

data = {key:{"mean":[], "std":[], "x":[]} for key in keys}

for mode in modes:
    path = ("../celmaWithBoussinesqApprox/BousCSDXMagFieldScanAr/"
            "visualizationPhysical/B0_{}/"
            "performance/performance{}.pickle").format(B0,mode)

    with open(path, "rb") as f:
        fig = pickle.load(f)

    axDown = fig.get_axes()[1]
    lines   = axDown.get_lines()
    for l in lines:
        label = l.get_label().replace("$","").\
                    replace(r"\mathrm{","").replace("}","").\
                    replace(r"\quad ","")

        yData = l.get_data()[1]

        # Ensure postivie numbers (bug in log writer)
        # Also makes all bars visible
        mean = yData.mean()
        if mean < 1:
            mean = 1
            std  = 0
        else:
            std = yData.std()

        data[label]["mean"].append(mean)
        data[label]["std"].append(std)

    plt.close(fig)

# Make segments
# +3 as we would like some space between the bars
lenKeys  = len(keys)+5
nBars    = len(modes)*lenKeys
xBarVals = tuple(range(nBars))

data["Arithmetics"       ]["x"] = xBarVals[0::lenKeys]
data["Communication"     ]["x"] = xBarVals[1::lenKeys]
data["Input/output"      ]["x"] = xBarVals[2::lenKeys]
data["Laplace inversions"]["x"] = xBarVals[3::lenKeys]
data["Time solver"       ]["x"] = xBarVals[4::lenKeys]

tickVals= data["Input/output"]["x"]

# Create the figure
fig, ax = plt.subplots(figsize = SizeMaker.standard(a=0.5, s=0.5))

for key in keys:
    d = data[key]
    ax.bar(d["x"], d["mean"], yerr=d["std"], label=key)

PlotHelper.makePlotPretty(ax)
ax.xaxis.grid(False)
ax.xaxis.set_ticks(tickVals)
ax.xaxis.set_ticklabels(("Initial\nphase", "Expand\nphase", "Linear\nphase", "Turbulent\nphase"))
ax.set_ylabel("$\%$")

# Move legend outside
handles, labels = ax.get_legend_handles_labels()

# Remove old legend
leg = ax.legend()
leg.remove()
fig.legend(handles,\
           labels ,\
           bbox_to_anchor=(1.05, 1.0),\
           loc="upper left",\
           borderaxespad=0.,\
           bbox_transform = ax.transAxes,\
           )

PlotHelper.savePlot(fig, "PercentagesBousScanB0{}".format(B0).replace(".","_")+".pdf")
