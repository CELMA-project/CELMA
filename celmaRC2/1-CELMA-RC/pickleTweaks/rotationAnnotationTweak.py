#!/usr/bin/env python

"""
Tweak of profile plots to be used in thesis.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import PlotHelper

picklePath = "../CSDXMagFieldScanAr/visualizationPhysical/B0_0.1/field2D/n-perpPol-2D-fluct-0.pickle"
with open(picklePath, "rb") as f:
    fig = pickle.load(f)

axes   = fig.get_axes()
perpAx = axes[0]
polAx  = axes[2]

color = None
arrowstyle = "simple"
txt = ""

# Circle arrow
coord = "data"
start = (- 0.06, 0.04 )
end   = (  0.00, 0.075)
connectionstyle = "arc3,rad=-0.3"
perpAx.annotate(txt,
                xytext=start, textcoords=coord,
                xy=end, xycoords=coord,
                arrowprops=dict(arrowstyle=arrowstyle,
                                color=color,
                                connectionstyle=connectionstyle,
                                ),
                )

# Straigth arrow
# NOTE: "axes fraction" equals ax.transAxes
coord = "axes fraction"
start = (0.85, 0.85)
end   = (0.5 , 0.85)
connectionstyle = "arc3,rad=0.0"
polAx.annotate(txt,
               xytext=start, textcoords=coord,
               xy=end, xycoords=coord,
               arrowprops=dict(arrowstyle=arrowstyle,
                               color=color,
                               connectionstyle=connectionstyle,
                               ),
               )

fileName = os.path.splitext(picklePath)[0] + "_rot.pdf"
PlotHelper.savePlot(fig, fileName)
