#!/usr/bin/env python

"""
Tweak of Boussinesq profile plots to be used in thesis.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np
from subprocess import Popen

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper, seqCMap3

folder = "../../B1-CELMA-RC/BousCSDXMagFieldScanAr/visualizationPhysical/B0_0.1/field1D"
for direction in ("radial", "parallel"):
    picklePath = os.path.join(folder,\
                              "mainFieldsBoussinesq-{}-1D-0.pickle".format(direction))
    with open(picklePath, "rb") as f:
        fig = pickle.load(f)

    lnAx, phiAx, nAx, omAx, jPaAx, nuiAx, uiAx, sAx, ueAx =\
            fig.get_axes()

    # Swap axes
    phiAx.set_position(omAx.get_position())
    omAx .set_position(nuiAx.get_position())

    fig.delaxes(lnAx)
    fig.delaxes(nuiAx)

    # Modify title position
    t = fig.texts[0]
    pos = list(t.get_position())
    pos[1] = 0.75
    t.set_position(pos)
    t.set_va("bottom")

    # Color adjust
    axes = (nAx, phiAx, jPaAx, omAx, uiAx, sAx, ueAx)
    colors = seqCMap3(np.linspace(0, 1, len(axes)))

    # Recolor the lines
    for ax, color in zip(axes, colors):
        line = ax.get_lines()[0]
        line.set_color(color)
        # NOTE: Using the designated setter gives a comparion error
        line._markerfacecolor = color
        line._markeredgecolor = color

        # To fix the legends, it seems like the easiest to do is to
        # recreate them
        handles, labels = ax.get_legend_handles_labels()
        leg = ax.legend()
        # Remake legend
        leg.remove()
        ax.legend(handles         ,\
                  labels          ,\
                  loc="upper left",\
                 )

    if direction == "radial":
        fileName = "B010RadBous.pdf"
    elif direction == "parallel":
        fileName = "B010ParBous.pdf"

    # Let pdfcrop do the cropping as "tigth" cuts some text
    PlotHelper.savePlot(fig, fileName, crop=False)
    Popen("pdfcrop {0} {0}".format(fileName), shell=True).wait()
