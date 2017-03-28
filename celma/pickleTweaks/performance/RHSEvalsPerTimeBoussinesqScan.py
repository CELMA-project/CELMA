#!/usr/bin/env python

"""
Plots the RHS evaluations per timestep for the Boussinesq B-scan.
"""

import pickle
import matplotlib.pylab as plt

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper

B0s     = (0.02, 0.04, 0.06, 0.08, 0.1)
modes   = ("Init", "Expand", "Linear", "Turbulence")

# Make place holders for the y-values
initMeans       = []
expandMeans     = []
linearMeans     = []
turbulenceMeans = []

initStds       = []
expandStds     = []
linearStds     = []
turbulenceStds = []

means = (initMeans, expandMeans, linearMeans, turbulenceMeans)
stds  = (initStds, expandStds, linearStds, turbulenceStds)

for B0 in B0s:
    for mode, mean, std in zip(modes, means, stds):
        path = ("../celmaWithBoussinesqApprox/BousCSDXMagFieldScanAr/"
                "visualizationPhysical/B0_{}/"
                "performance/performance{}.pickle").format(B0,mode)

        try:
            with open(path, "rb") as f:
                fig = pickle.load(f)

            axUp = fig.get_axes()[0]
            _, RHSEvals = axUp.get_lines()[0].get_data()

            mean.append(RHSEvals.mean())
            std .append(RHSEvals.std() )

            plt.close(fig)
        except FileNotFoundError:
            mean.append(0)
            std .append(0)

# Make segments
# +1 as we would like some space between the bars
lenMode         = len(modes)+1
nBars           = len(B0s)*lenMode
xBarVals        = tuple(range(nBars))

initXvals       = xBarVals[0::lenMode]
expandXvals     = xBarVals[1::lenMode]
linearXvals     = xBarVals[2::lenMode]
turbulenceXvals = xBarVals[3::lenMode]
tickVals        = tuple(val - 0.5 for val in linearXvals)

# Create the figure
fig, ax = plt.subplots(figsize = SizeMaker.standard(a=0.5, s=0.5))

ax.bar(initXvals, initMeans,\
        yerr=initStds, label="Initial phase")
ax.bar(expandXvals, expandMeans,\
        yerr=expandStds, label="Expand phase")
ax.bar(linearXvals, linearMeans,\
        yerr=linearStds, label="Linear phase")
ax.bar(turbulenceXvals, turbulenceMeans,\
        yerr=turbulenceStds, label="Turbulent phase")

PlotHelper.makePlotPretty(ax)
ax.xaxis.grid(False)
ax.xaxis.set_ticks(tickVals)
ax.xaxis.set_ticklabels(B0s)
ax.set_xlabel("$B_0 [T]$")
ax.set_ylabel("RHS iterations\nper time step")

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

PlotHelper.savePlot(fig, "RHSEvalsPerTimeBoussinesqScan.pdf")
