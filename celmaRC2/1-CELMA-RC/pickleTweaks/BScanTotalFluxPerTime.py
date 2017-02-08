#!/usr/bin/env python

"""
Gives blobs per unit time for the B-Scan.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper

sigMul = 3
scans = ("B0_0.06", "B0_0.08", "B0_0.1")
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
    path = "../CSDXMagFieldScanAr/visualizationPhysical/{}/totalFluxes/".format(scan)

    # Obtain the time
    fileName = os.path.join(path, "totalFluxes.pickle")
    with open(fileName, "rb") as f:
        fig = pickle.load(f)

    parallelAx, perpAx = fig.get_axes()

    time, _ = perpAx.get_lines()[0].get_data()
    lastT = time[-1]
    curScan = float(scan[4:])

    # Get the numbers from the first axis
    t = parallelAx.texts[0].get_text()
    totParEl  = getNumber(t, ind = 1)
    totParIon = getNumber(t, ind = 2)

    t = perpAx.texts[0].get_text()
    totPerp = getNumber(t, ind = 1)

    elPerTime  .append(totParEl/lastT)
    ionsPerTime.append(totParIon/lastT)
    perpPerTime.append(totPerp/lastT)

    scanVal.append(curScan)


fig, (parAx, perpAx) = plt.subplots(ncols = 2, figsize = SizeMaker.array(2,1))
parAx.plot(scanVal, elPerTime, ls="", alpha = 0.7,\
           marker="o", ms=10, label="$\mathrm{Electrons}$")
parAx.plot(scanVal, ionsPerTime, ls="", alpha = 0.7,\
           marker="d", ms=10, label="$\mathrm{Ions}$")
perpAx.plot(scanVal, perpPerTime, "g", ls="", alpha = 0.7,\
            marker="s", ms=10, label="$\mathrm{Perpendicular}$")
parAx .set_ylabel(r"$[\mathrm{s}^{-1}]$")
perpAx.set_ylabel(r"$[\mathrm{s}^{-1}]$")
parAx .set_xlabel(r"$B [\mathrm{T}]$")
perpAx.set_xlabel(r"$B [\mathrm{T}]$")

parAx .set_title(r"$\mathrm{Parallel}$")
perpAx.set_title(r"$\mathrm{Perpendicular}$")

fig.suptitle(r"$\mathrm{Average\;particle\;flux\;per\;second}$", y=1.1)
PlotHelper.makePlotPretty(parAx , rotation = 45)
PlotHelper.makePlotPretty(perpAx, rotation = 45)

fig.subplots_adjust(wspace=0.75)
PlotHelper.savePlot(fig, "BScanTotalFlux.pdf")
