#!/usr/bin/env python

"""
Gives blobs per unit time for the nn-Scan.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper, plotNumberFormatter

scans  = (\
          "nn_9.9e+20"              ,\
          "nn_4e+19"                ,\
          "nn_1.5e+19"              ,\
          "nn_2.5e+18"              ,\
          "nn_6.666666666666668e+18",\
         )

elPerTime   = []
ionsPerTime = []
perpPerTime = []
scanVal     = []

#{{{getNumber
def getNumber(t, ind):
    """
    Retruns the string in the total fluxes.

    Parameters
    ----------
    t : str
        String to extract the number from
    ind : int
        Index to investigate after the first split.
    """

    return eval(t.split("=")[ind].split(";")[0].replace("$","") .\
                                                replace("\\","").\
                                                replace("{","") .\
                                                replace("}","") .\
                                                replace("cdot 10^","e").
                                                replace("cdot","*"))
#}}}

for scan in scans:
    path = "../../CSDXNeutralScanAr/visualizationPhysical/{}/totalFluxes/".format(scan)

    # Obtain the time
    fileName = os.path.join(path, "totalFluxes.pickle")
    with open(fileName, "rb") as f:
        fig = pickle.load(f)

    parallelAx, perpAx = fig.get_axes()

    time, _ = perpAx.get_lines()[0].get_data()
    lastT = time[-1]
    curScan = float(scan[3:])

    # Get the numbers from the first axis
    t = parallelAx.texts[0].get_text()
    totParEl  = getNumber(t, ind = 1)
    totParIon = getNumber(t, ind = 2)

    t = perpAx.texts[0].get_text()
    totPerp = getNumber(t, ind = 1)

    elPerTime  .append(totParEl/lastT)
    ionsPerTime.append(totParIon/lastT)
    perpPerTime.append(totPerp/lastT)

    scanVal.append(plotNumberFormatter(curScan, None))

plt.close("all")

# Perp values
fig, (perpAx, parAx) = plt.subplots(ncols=2, figsize = (5,1))
yTicks   = tuple(r"${}\; \mathrm{{T}}$".format(B) for B in scanVal)
xBarVals = tuple(range(len(perpPerTime)))

# Find the lowest value
lowestVal = np.min(perpPerTime)
# Pad with 5%
lowestVal = 0.95*lowestVal

# Get the color cycler
# http://stackoverflow.com/questions/13831549/get-matplotlib-color-cycle-state
prop_cycler = perpAx._get_lines.prop_cycler
# Plot
perpRects = perpAx.barh(xBarVals                        ,\
                        perpPerTime                     ,\
                        left = lowestVal                ,\
                        color=next(prop_cycler)["color"],\
                       )

# Set text
for nr, (rect, txt) in enumerate(zip(perpRects, yTicks)):
    width  = rect.get_width()
    height = rect.get_y()
    # Add 4% padding
    x      = width*1.04 + lowestVal
    y      = height + rect.get_height()/2
    perpAx.text(x, y, txt, ha="left", va="center")

perpAx.set_xlabel(r"$[\mathrm{s}^{-1}]$")
perpAx.set_title(r"$\mathrm{Perpendicular}$")

perpAx.spines["top"]  .set_visible(False)
perpAx.spines["left"] .set_visible(False)
perpAx.spines["right"].set_visible(False)
perpAx.get_yaxis().set_visible(False)

PlotHelper.makePlotPretty(perpAx, xbins=5, legend=False, rotation=45)
perpAx.grid(False)


# Par values
xBarVals = tuple(range(len(perpPerTime)*2))
ionXBarVals = xBarVals[1::2]
elXBarVals  = xBarVals[0::2]

# Find the lowest value
lowestVal = np.min(perpPerTime)
# Pad with 5%
lowestVal = 0.95*lowestVal

# Plot
elRects = parAx.barh(elXBarVals                      ,\
                     elPerTime                       ,\
                     left = lowestVal                ,\
                     color=next(prop_cycler)["color"],\
                    )
ionRects = parAx.barh(ionXBarVals                     ,\
                      ionsPerTime                     ,\
                      left = lowestVal                ,\
                      color=next(prop_cycler)["color"],\
                     )

# Set text
for nr, (rect, txt) in enumerate(zip(ionRects, yTicks)):
    width  = rect.get_width()
    height = rect.get_y()
    # Add 4% padding
    x      = width*1.04 + lowestVal
    y      = height
    parAx.text(x, y, txt, ha="left", va="center")

parAx.set_xlabel(r"$[\mathrm{s}^{-1}]$")
parAx.set_title(r"$\mathrm{Parallel}$")

parAx.spines["top"]  .set_visible(False)
parAx.spines["left"] .set_visible(False)
parAx.spines["right"].set_visible(False)
parAx.get_yaxis().set_visible(False)

PlotHelper.makePlotPretty(parAx, xbins=5, legend=False, rotation=45)
parAx.grid(False)

# Make legend manually
handles = (perpRects[0], ionRects[0], elRects[0])
labels  = ("$\mathrm{Both}$"     ,\
           "$\mathrm{Ions}$"     ,\
           "$\mathrm{Electrons}$",\
           )

fig.legend(handles                         ,\
           labels                          ,\
           bbox_to_anchor=(1.25, 1.0)      ,\
           loc="upper left"                ,\
           borderaxespad=0.                ,\
           bbox_transform = parAx.transAxes,\
           )

fig.suptitle(r"$\mathrm{Average\;particle\;flux\;per\;second}$", y=1.3)

PlotHelper.savePlot(fig, "nnScanTotalFlux.pdf")
