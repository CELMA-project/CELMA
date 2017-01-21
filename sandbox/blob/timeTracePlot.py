#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 13:48:38 2017

@author: mmag
"""

import pickle
import matplotlib.pylab as plt

val = "0.006Flux"
with open("{}.pickle".format(val), "rb") as f:
    rf = pickle.load(f)

keys = sorted(rf.keys())

for key in keys:
    fig, ax = plt.subplots()
    time = rf[key]["time"]
    lastKey = "n"
    if "Flux" in val:
        lastKey += "RadialFlux"
    var  = rf[key][lastKey]
    std  = var.std()
    ax.plot(time, var)
    # One std
    ax.plot((time[0],time[-1]), (std,std))
    ax.plot((time[0],time[-1]), (2*std,2*std))
    ax.plot((time[0],time[-1]), (3*std,3*std))
    rho = "{:.3f}".format(float(key.split(",")[0]))
    ax.set_title(rho)
    ax.grid()

    # Max windows
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

    fig.savefig("{}_{}.png".format(val,rho))

plt.show()
