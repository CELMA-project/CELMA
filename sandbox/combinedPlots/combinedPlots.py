#!/usr/bin/env python

"""Post-processor test for timeTraces"""

import pickle
import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.timeTrace import PlotTimeTrace
from CELMAPy.unitsConverter import UnitsConverter
import numpy as np


from matplotlib.gridspec import GridSpec
import matplotlib.pylab as plt

plt.rc("text", usetex=False)



# NOTE: time slicing happens in the collect routine
with open("./timeTraces.pickle", "rb") as f:
    tt = pickle.load(f)

# Mess with the signal
size = 5000
for nr, key in enumerate(tt.keys()):
    tt[key]["n"]    = (2*np.random.random((size,))-1)/(nr+1)
    tt[key]["time"] = np.array(range(size))

# Oh sh..forgot uc object
# Fix here
path ="CSDXMagFieldScan/nout_1000_timestep_1/geom_Lx_19.8078_geom_Ly_693.2731_input_B0_0.04_switch_forceAddNoise_False_switch_includeNoise_False_switch_saveTerms_False_switch_useHyperViscAzVortD_True_tag_CSDXMagFieldScan-2-linearPhase1_0/"
uc =UnitsConverter(path)
# DONE

savePath = "."
mode = "fluct"

plotSuperKwargs = {\
                    "showPlot"     : True,\
                    "savePlot"     : False,\
                    "savePath"     : savePath,\
                    "savePathFunc" : None,\
                    "extension"    : None,\
                    "dmp_folders"  : None,\
                    "sliced"       : False,\
                   }



pltSize = (20,20)
fig = plt.figure(figsize=pltSize)

gs = GridSpec(len(tt.keys()), 3)
axes = []
# Make first ax
firstAx = fig.add_subplot(gs[0])
firstAx.grid(True)
axes.append(firstAx)
# FIXME: Don't know where to share xax, maybe don't need
# Make the rest of the axes
#NOTE: The gridspec works like a book:
#       0, 1, 2
#       3, 4, 5
# In other words, time traces will be plotted on 0,3,6. ie n%3=0
# PDF will be plotted on n%3=1
# PSD will be plotted on n%3=2

#{{{timeTrace
# Slice
tSlice = slice(10,50)
timeTraceVar  = []
timeTraceTime = []
for key in tt.keys():
    timeTraceVar .append(tt[key]["n"]   [tSlice])
    timeTraceTime.append(tt[key]["time"][tSlice])

timeTraceVar  = tuple(timeTraceVar)
timeTraceTime = tuple(timeTraceTime)
#}}}

YOU ARE HERE
#{{{PDF
# Slice
tSlice = slice(10,50)
timeTraceVar  = []
timeTraceTime = []
for key in tt.keys():
    timeTraceVar .append(tt[key]["n"]   [tSlice])
    timeTraceTime.append(tt[key]["time"][tSlice])

timeTraceVar  = tuple(timeTraceVar)
timeTraceTime = tuple(timeTraceTime)
#}}}

# Get max min of all
timeTraceMin = np.min(timeTraceVar)
timeTraceMax = np.max(timeTraceVar)

for nr in range(9):
    ax = fig.add_subplot(gs[nr])
    ax.grid(True)
    axes.append(ax)
    if nr%3==0:
        ax.plot(timeTraceTime[nr//3], timeTraceVar[nr//3])
        ax.set_ylim((timeTraceMin, timeTraceMax))
    elif nr%3==1:
        # PDF
        ax.plot([nr])
    elif nr%3==2:
        #PSD
        ax.plot([nr])

plt.show()

# # Plot
# # FIXME: Notice own uc
# ptt = PlotTimeTrace(uc         ,\
#                     **plotSuperKwargs)
# ptt.setData(tt, mode, timeAx=False)
# ptt.plotSaveShowTimeTrace()
