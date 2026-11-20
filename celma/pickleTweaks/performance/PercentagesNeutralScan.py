#!/usr/bin/env python

"""
Plots the percentages for the neutral-scan.
"""

import os
import pickle
import sys
from typing import TypedDict


class PerformanceData(TypedDict):
    """Collected means, deviations, and bar positions."""

    mean: list[float]
    std: list[float]
    x: list[int] | tuple[int, ...]


import matplotlib.pylab as plt

# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper, SizeMaker

B0 = 0.06
nn = 1.5e19
modes = ("Init", "Expand", "Linear", "Turbulence")

# Make place holders for the data
keys = (
    "Arithmetics",
    "Communication",
    "Input/output",
    "Laplace inversions",
    "Time solver",
)

data: dict[str, PerformanceData] = {
    key: {"mean": [], "std": [], "x": []} for key in keys
}

for mode in modes:
    path = (
        "../../CSDXNeutralScanAr/"
        "visualizationPhysical/nn_{}/"
        "performance/performance{}.pickle"
    ).format(nn, mode)

    with open(path, "rb") as f:
        fig = pickle.load(f)

    axDown = fig.get_axes()[1]
    lines = axDown.get_lines()
    for l in lines:
        label = (
            l.get_label()
            .replace("$", "")
            .replace(r"\mathrm{", "")
            .replace("}", "")
            .replace(r"\quad ", "")
        )

        yData = l.get_data()[1]

        # Ensure postivie numbers (bug in log writer)
        # Also makes all bars visible
        mean = yData.mean()
        if mean < 1:
            mean = 1
            std = 0
        else:
            std = yData.std()

        data[label]["mean"].append(mean)
        data[label]["std"].append(std)

    plt.close(fig)

# Make segments
# +3 as we would like some space between the bars
lenKeys = len(keys) + 5
nBars = len(modes) * lenKeys
xBarVals = tuple(range(nBars))

data["Arithmetics"]["x"] = xBarVals[0::lenKeys]
data["Communication"]["x"] = xBarVals[1::lenKeys]
data["Input/output"]["x"] = xBarVals[2::lenKeys]
data["Laplace inversions"]["x"] = xBarVals[3::lenKeys]
data["Time solver"]["x"] = xBarVals[4::lenKeys]

tickVals = data["Input/output"]["x"]

# Create the figure
fig, ax = plt.subplots(figsize=SizeMaker.standard(a=0.5, s=0.5))

for key in keys:
    d = data[key]
    ax.bar(d["x"], d["mean"], yerr=d["std"], label=key)

PlotHelper.makePlotPretty(ax)
ax.xaxis.grid(False)
ax.xaxis.set_ticks(tickVals)
ax.xaxis.set_ticklabels(
    ("Initial\nphase", "Expand\nphase", "Linear\nphase", "Turbulent\nphase")
)
ax.set_ylabel(r"$\%$")

# Move legend outside
handles, labels = ax.get_legend_handles_labels()

# Remove old legend
leg = ax.legend()
leg.remove()
fig.legend(
    handles,
    labels,
    bbox_to_anchor=(1.05, 1.0),
    loc="upper left",
    borderaxespad=0.0,
    bbox_transform=ax.transAxes,
)

PlotHelper.savePlot(
    fig, "PercentagesNeutralScanB0{}Nn{}".format(B0, nn).replace(".", "_") + ".pdf"
)
