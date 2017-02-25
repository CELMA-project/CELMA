#!/usr/bin/env python

"""
Tweak profile comparison between non-Boussinesq and Boussinesq.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np
from subprocess import Popen

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper, seqCMap3

folder = "."

for direction in ("Rad", "Par"):
    # Store in dict in order to avoid arbitrary order
    axDict = {}

    # Open the non-Boussinesq plot
    picklePath = os.path.join(folder,\
                              "B010{}.pickle".format(direction))
    with open(picklePath, "rb") as f:
        fig = pickle.load(f)


    axes = fig.get_axes()

    for ax in axes:
        line = ax.get_lines()[0]
        line.set_alpha(0.7)
        line.set_markevery(7)
        line.set_markersize(5)

        handles, labels = ax.get_legend_handles_labels()
        # Store into dict
        ind = labels[0][:5]
        axDict[ind] = {}
        axDict[ind]["ax"]     = ax
        axDict[ind]["handle"] = handles[0]
        axDict[ind]["label"]  = labels[0].split(" ")[0] + r" $\mathrm{F}$"
        axDict[ind]["color"]  = line.get_color()

        # Remove old legend
        leg = ax.legend()
        leg.remove()

    # Open the Boussinesq plot
    picklePath = os.path.join(folder,\
                              "B010{}Bous.pickle".format(direction))
    with open(picklePath, "rb") as f:
        bFig = pickle.load(f)

    bAxes = bFig.get_axes()

    for bAx in bAxes:
        line = bAx.get_lines()[0]
        x, y = line.get_data()

        # Get the legend label
        _, bl = bAx.get_legend_handles_labels()

        # Get the index
        ind = bl[0][:5]

        # Fix the legend
        bl = bl[0].split(" ")[0] + r" $\mathrm{B}$"


        axDict[ind]["ax"].plot(x, y,\
                               markevery = 7, markersize = 5, marker = "^",\
                               ls = ":",\
                               color = axDict[ind]["color"], alpha = 0.7,\
                               label = bl)

        # Update the handles and labels, and recreate them
        handles, labels = axDict[ind]["ax"].get_legend_handles_labels()

        # Remake legend
        leg = axDict[ind]["ax"].legend()
        leg.remove()
        axDict[ind]["ax"].legend((axDict[ind]["handle"], handles[1]),\
                                 (axDict[ind]["label"] , labels [1]),\
                                 loc="upper left",\
                                )

        # Autoscale
        axDict[ind]["ax"].autoscale(enable=True, axis="y", tight=True)

    if direction == "Rad":
        fileName = "B010RadBousCompare.pdf"
    elif direction == "Par":
        fileName = "B010ParBousCompare.pdf"

    # Let pdfcrop do the cropping as "tigth" cuts some text
    PlotHelper.savePlot(fig, fileName, crop=False)
    Popen("pdfcrop {0} {0}".format(fileName), shell=True).wait()
