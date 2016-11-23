#!/usr/bin/env python

"""
Scripts which plots the valid region of drift approximation.

NOTE: Although we find the region by inversion, the zero we find is
      proporitonal to B/m_i.
      To convert to another gas at another B, multiply the y-axis with
      (m_iNew/m_iOld)/(BNew/BOld)
"""

from scipy.optimize import fsolve
import scipy.constants as cst
import numpy as np
import matplotlib.pylab as plt

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../celma/common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotHelpers import plotNumberFormatter, PlotHelper

showPlot  = True
extention = "pdf"
# Can be Ar or H
gas = "H"
B   = 1e-1

me = cst.m_e
q  = cst.e
e0 = cst.epsilon_0
pi = np.pi

if gas == "H":
    mi = cst.m_p
elif gas == "Ar":
    mi = 39.948*cst.u

TeEV = np.linspace(1, 15, 20)
TeJ  = TeEV*cst.eV

def nuEIOverOmCIMinusOne(n, Te, B, mi):
    """
    Calculates $\nu_ei/\omega_{ci} = 1$

    We have that
    >>> coulombLog = log(12.0*pi*(e0*Te0)**(3/2)/(n0*q**3))
    >>> vThE = (2.0*Te0/me)**0.5
    >>> nuEI = mi*ne*q**3*coulombLog/(4*B*e0**2*me**2*pi*vThE**3.0)
    >>> omCI = B*q/mi

    Parameters
    ----------
    n : float
        The density
    Te : float
        The electron temperature in J
    B : float
        The magnetic field
    mi : float
        The ion mass
    """
    shouldBeZero =\
            (mi*n*(q**3.0)*np.log(12.0*pi*((e0*Te)**(3.0/2.0))/(np.sqrt(n)*q**3.0)))/\
            (4.0*B*(e0**2.0)*(me**2.0)*pi*((2.0*Te/me)**(3.0/2.0)))\
            - 1
    return shouldBeZero

fig, ax = plt.subplots()

n = np.zeros(TeJ.shape)
for i in range(len(TeJ)):
    # Calculates for what n $\nu_ei/\omega_{ci} = 1$
    n[i] = fsolve(nuEIOverOmCIMinusOne, x0=1e16, args=(TeJ[i], B, mi))

n = tuple(n)

ax.fill_between(TeEV, n, facecolor="b", edgecolor="b", alpha=0.5)
ax.set_xlim((np.min(TeEV), np.max(TeEV)))
ax.set_xlabel("Te [eV]")
ax.set_ylabel("n [$m^{-3}$]")
fig.suptitle("{} at B$={}$".format(gas, B))

xAnnotate = np.max(TeEV)/3
yAnnotate = np.max(n)*3/4

size = 30
ha   = "center"

ax.annotate("Drift\napproximation\nbreaks",\
            # Location of annotation
            xy=(xAnnotate, yAnnotate),\
            # Location of text
            xycoords="data",\
            xytext=(xAnnotate, yAnnotate),\
            textcoords="data",
            size=size,\
            va="center",\
            ha=ha,\
            bbox=dict(boxstyle="round", fc="w", alpha=0.5),\
            # An arrow is not used here
            )

xAnnotate = np.max(TeEV)*3/4
yAnnotate = np.max(n)/4

ax.annotate("Drift\napproximation\nvalid",\
            # Location of annotation
            xy=(xAnnotate, yAnnotate),\
            # Location of text
            xycoords="data",\
            xytext=(xAnnotate, yAnnotate),\
            textcoords="data",
            size=size,\
            va="center",\
            ha=ha,\
            bbox=dict(boxstyle="round", fc="w", alpha=0.5),\
            # An arrow is not used here
            )

PlotHelper.makePlotPretty(ax, rotation=45)

if extention:
    PlotHelper.savePlot(fig, "validRegion.{}".format(extention))

if showPlot:
    plt.show()
