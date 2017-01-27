#!/usr/bin/env python

"""
Example on how to edit a pickled plot.
"""

import pickle
import matplotlib.pylab as plt

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper

path = "../CSDXMagFieldScanAr/visualizationPhysical/all/growthRates/growthRatesAnalytic-rho-00388-z-0.7-B0.pickle"
with open(path, "rb") as f:
    fig = pickle.load(f)

axUp = fig.get_axes()[0]

# Set the alpha
lines = axUp.get_lines()

for line in lines:
    line.set_alpha(0.7)
    line.set_markersize(7)

handles, labels = axUp.get_legend_handles_labels()

leg = axUp.legend()

# Remove old legend
leg.remove()

fig.legend(handles,\
           labels,\
           bbox_to_anchor=(1.05, 1.0),\
           loc="upper left",\
           borderaxespad=0.,\
           bbox_transform = axUp.transAxes,\
           )

# Tweak the size a bit
fig.set_figwidth(2.0)
fig.set_figheight(4.0)

# Modify title position
t = fig.texts[0]
pos = list(t.get_position())
pos[1] = 1.05
t.set_transform(axUp.transAxes)
t.set_position(pos)
t.set_va("bottom")

PlotHelper.savePlot(fig, "test2.pdf")
