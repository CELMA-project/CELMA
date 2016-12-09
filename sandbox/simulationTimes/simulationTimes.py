#!/usr/bin/env python

"""
Plots the simulation times
"""

import numpy as np
from numpy import log, sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
import scipy.constants as cst
import datetime

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../celma/common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotHelpers import plotNumberFormatter, PlotHelper

showPlots = True
savePlots = "png"

# Dict building on celmaRC/CSDXMagFieldScanAr/conclusion.txt
CSDXArTimes =\
{
    "Initial phase":\
    {
     0.02 : {"Time":"04h43m", "steps":4000, "Success":True},\
     0.04 : {"Time":"15h30m", "steps":4000, "Success":True},\
     0.06 : {"Time":"19h41m", "steps":4000, "Success":True},\
     0.08 : {"Time":"26h21m", "steps":4000, "Success":True},\
     0.1  : {"Time":"26h51m", "steps":4000, "Success":True},\
    },\
    "Expand phase":\
    {
     0.02 : {"Time":"10h21m", "steps":100, "Success":True },\
     0.04 : {"Time":"23h43m", "steps":100, "Success":True },\
     0.06 : {"Time":"72h00m", "steps": 50, "Success":False},\
     0.08 : {"Time":"20h24m", "steps": 50, "Success":False},\
     0.1  : {"Time":"22h14m", "steps": 50, "Success":False},\
    },\
    "Linear phase":\
    {
     0.02 : {"Time":"32h26m", "steps":1000     , "Success":True },\
     0.04 : {"Time":"72h00m", "steps":4881-4100, "Success":False},\
     0.06 : {"Time":"72h00m", "steps":4917-4050, "Success":False},\
     0.08 : {"Time":"00h00m", "steps":1000     , "Success":False},\
     0.1  : {"Time":"00h00m", "steps":1000     , "Success":False},\
    },\
}

CSDXArBTimes =\
{
    "Initial phase":\
    {
     0.02 : {"Time":"00h15m", "steps":6000, "Success":False},\
     0.04 : {"Time":"00h14m", "steps":4000, "Success":True},\
     0.06 : {"Time":"03h53m", "steps":4000, "Success":True},\
     0.08 : {"Time":"07h17m", "steps":4000, "Success":True},\
     0.1  : {"Time":"12h47m", "steps":4000, "Success":True},\
    },\
    "Expand phase":\
    {
     0.02 : {"Time":"00h16m", "steps":100, "Success":True},\
     0.04 : {"Time":"00h14m", "steps":100, "Success":True},\
     0.06 : {"Time":"00h05m", "steps":100, "Success":True},\
     0.08 : {"Time":"09h06m", "steps":100, "Success":True},\
     0.1  : {"Time":"15h32m", "steps":100, "Success":True},\
    },\
    "Linear phase":\
    {
     0.02 : {"Time":"00h00m", "steps":1000     , "Success":False},\
     0.04 : {"Time":"04h46m", "steps":1000     , "Success":True},\
     0.06 : {"Time":"72h00m", "steps":4942-4100, "Success":False},\
     0.08 : {"Time":"72h00m", "steps":4775-4100, "Success":False},\
     0.1  : {"Time":"72h00m", "steps":4381-4100, "Success":False},\
    },\
    "Turbulence phase":\
    {
     0.02 : {"Time":"39h43m", "steps":5000, "Success":True },\
     0.04 : {"Time":"00h00m", "steps":5000, "Success":False},\
     0.06 : {"Time":"00h00m", "steps":5000, "Success":False},\
     0.08 : {"Time":"00h00m", "steps":5000, "Success":False},\
     0.1  : {"Time":"00h00m", "steps":5000, "Success":False},\
    },\
}

colors  = ["c", "m", "k", "r", "g", "b"]

def plotTimes(times, title, scanParameter, includeTotal=True):
    """
    Plot the runtimes

    Parameters
    ----------
    times : dict
        Dictionary of all the times
    title : str
        Title of the plot
    scanParameter: str
        String of the scan parameter
    includeTotal : bool
        Whether the summed time should be used or not
    """
    fig, ax = plt.subplots()
    legends = {"scat":[], "legend":[]}
    count   = 0

    if includeTotal:
        total = {}

    for runType, runTypeDict in times.items():
        x = list()
        y = list()
        for scan, scanDict in runTypeDict.items():
            x.append(scan)
            seconds = datetime.timedelta(hours   = int(scanDict["Time"][0:2]),\
                                         minutes = int(scanDict["Time"][3:5]))\
                               .total_seconds()

            steps          = scanDict["steps"]
            minutesPrSteps = (seconds/60)/steps
            y.append(minutesPrSteps)

            if includeTotal:
                if not(scan in total.keys()):
                    # NOTE: We would like sum(seconds)/sum(steps)
                    #       rather than sum(seconds/steps)
                    total[scan] = {"seconds":0, "steps":0}
                total[scan]["seconds"] += seconds
                total[scan]["steps"]   += steps

        x = tuple(x)
        y = tuple(y)
        scat = ax.scatter(x, y, color=colors[count], s=100, alpha=0.7)
        legends["legend"].append(runType)
        legends["scat"  ].append(scat)
        count += 1

    if includeTotal:
        minutesPrStepsTotal = {}
        for key in total.keys():
            minutesPrStepsTotal[key] =\
               (total[key]["seconds"]/60)/total[key]["steps"]

        x, y = zip(*minutesPrStepsTotal.items())
        scat = ax.scatter(x, y,\
                          color=colors[count], s=100, alpha=0.7)
        legends["legend"].append("Total time/total steps")
        legends["scat"  ].append(scat)


    legends["legend"] = tuple(legends["legend"])
    legends["scat"  ] = tuple(legends["scat"  ])


    ax.set_xlabel(scanParameter)
    ax.set_ylabel("Minutes/iteration")

    leg = ax.legend(legends["scat"  ],\
                    legends["legend"],\
                    scatterpoints = 1,\
                    loc = "best"     ,\
                    fancybox = True
                   )

    leg.get_frame().set_alpha(0.5)

    ax.grid()

    fig.suptitle(title)

    plt.tight_layout()

if __name__ == "__main__":
    plotTimes(CSDXArTimes, "CSDX simulations with Ar", "$B_0$")
    plotTimes(CSDXArBTimes, "CSDX simulations with Ar using Boussinesq", "$B_0$")

    plt.show()
