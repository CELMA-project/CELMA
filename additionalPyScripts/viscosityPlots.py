#!/usr/bin/env python

"""
Scripts which plots the analytic viscosity coefficients.

These are given in for example
Helander, P. and Sigmar, D.J.
Collisional Transport in Magnetized Plasmas
"""

from parameters import Parameters
import numpy as np
from numpy import log, sqrt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.constants as cst

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../celmaRC2/common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import SizeMaker, PlotHelper, plotNumberFormatter

showPlots = False
savePlots = "pdf"

#{{{Common for simple plots
# Set the plotting style
title_height = 1.02
fig_no = 1
plt_size_n  = SizeMaker.standard(w=2.7, a=1)
plt_size_t = SizeMaker.standard(w=2.3, a=0.9)
title_size = 30
#}}}}

#{{{setNScanAxes
def setNScanAxes(fig, ax, yLabel, Te0EV):
    ax.set_xlabel(r"$n [\mathrm{m}^{-3}]$")
    ax.set_ylabel(yLabel)
    title = r"$Te=${} $[\mathrm{{eV}}]$".format(plotNumberFormatter(Te0EV, None))
    fig.suptitle(title, y=title_height)
    PlotHelper.makePlotPretty(ax, rotation=90, ybins=4)

    # Remove old legend
    leg = ax.legend()
    leg.remove()
    # Put the legend outside
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles,\
               labels,\
               bbox_to_anchor=(1.05, 1.0),\
               loc="upper left",\
               borderaxespad=0.,\
               bbox_transform = ax.transAxes,\
               )

    # Includes the xlabel if outside
    fig.tight_layout(rect=[0.001,0.001,1,1])
#}}}

#{{{setTScanAxes
def setTScanAxes(fig, ax, yLabel, n0):
    # Set axis label
    ax.set_xlabel(r"$Te [\mathrm{eV}]$")
    ax.set_ylabel(yLabel)
    title = r"$n=${} $[\mathrm{{m}}^{{-3}}]$".format(plotNumberFormatter(n0, None))
    fig.suptitle(title, y=title_height)
    PlotHelper.makePlotPretty(ax, rotation=90, ybins=4)

    # Remove old legend
    leg = ax.legend()
    leg.remove()
    # Put the legend outside
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles,\
               labels,\
               bbox_to_anchor=(1.05, 1.0),\
               loc="upper left",\
               borderaxespad=0.,\
               bbox_transform = ax.transAxes,\
               )

    # Includes the xlabel if outside
    fig.tight_layout(rect=[0.001,0.001,1,1])
#}}}

mi = cst.m_p
p = Parameters(mi)

B0s   = [1, 1e-1, 1e-2]
Ti0EV = 7.0
n0s    = np.linspace(1e16, 1e18, 200)
Te0sEV = np.linspace(7, 21, 200)

B0LineType = ["solid", "dashed", "dotted", "dash_dot"]
# One color for each eta
colors = cm.viridis(np.linspace(0, 1, 5))

# Make axes
figNuEINScan = plt.figure(0, figsize = plt_size_n)
figNuEITScan = plt.figure(1, figsize = plt_size_t)
figEtaINScan = plt.figure(2, figsize = plt_size_n)
figEtaITScan = plt.figure(3, figsize = plt_size_t)
figEtaENScan = plt.figure(4, figsize = plt_size_n)
figEtaETScan = plt.figure(5, figsize = plt_size_t)
axNuEINScan  = figNuEINScan.add_subplot(111)
axNuEITScan  = figNuEITScan.add_subplot(111)
axEtaENScan  = figEtaENScan.add_subplot(111)
axEtaINScan  = figEtaINScan.add_subplot(111)
axEtaETScan  = figEtaETScan.add_subplot(111)
axEtaITScan  = figEtaITScan.add_subplot(111)

etaE = {}

for B0, ls in zip(B0s, B0LineType):
    omCI = p.calcOmCI(B0)
    omCE = p.calcOmCE(B0)
    N    = 3.0

    # Convert to J
    Ti0 = Ti0EV*cst.eV

    #{{{n scan
    # Fixed Te0
    Te0EV = 14.0
    # Convert to J
    Te0 = Te0EV*cst.eV

    # Calculate normalization variables
    cS   = p.calcCS(Te0, Ti0, N)
    rhoS = p.calcRhoS(cS, omCI)

    # Caclulate collisionalities
    coulombLog = p.calcCoulombLog(Te0, n0s)
    vThE = p.calcVThE(Te0)
    vThI = p.calcVThI(Ti0)
    nuEI = p.calcNuEI(n0s, coulombLog, vThE)
    nuII = p.calcNuII(n0s, coulombLog, vThI)
    tauI = 1.0/(sqrt(2.0)*nuII)
    tauE = 1.0/nuEI

    # Caclulate ion viscosities
    etaI = p.calcEtaI(n0s, Ti0, tauI, omCI)
    sortedEtaIKeys = list(etaI.keys())
    sortedEtaIKeys.sort()

    # Caclulate electron viscosities
    etaE = p.calcEtaE(n0s, Te0, tauE, omCE)
    sortedEtaEKeys = list(etaE.keys())
    sortedEtaEKeys.sort()

    # Normalization
    nuEI /= omCI
    nuII /= omCI

    for key in sortedEtaIKeys:
        etaI[key] /= (mi*n0s*rhoS*cS)

    for key in sortedEtaEKeys:
        etaE[key] /= (mi*n0s*rhoS*cS)

    #{{{nuEI plot
    axNuEINScan.plot(n0s,\
                     nuEI,\
                     color = colors[0],\
                     ls    = ls,\
                     label = "$B_0={}$".format(B0))

    # Set axis
    setNScanAxes(figNuEINScan, axNuEINScan, r"$\nu_{ei}/\omega_{ci}$", Te0EV)

    if savePlots:
        PlotHelper.savePlot(figNuEINScan, "niEINScan.{}".format(savePlots))
    #}}}

    #{{{etaI plot
    for nr, key in enumerate(sortedEtaIKeys):
        if  B0 == B0s[0]:
            label = "${}$".format(key)
        else:
            label = None
        axEtaINScan.plot(n0s,\
                         etaI[key],\
                         color = colors[nr],\
                         ls    = ls,\
                         label = label)

    # Set axis
    axEtaINScan.set_yscale("log")
    setNScanAxes(figEtaINScan, axEtaINScan, r"$\eta_i/m_i n \rho_s c_s$", Te0EV)

    if savePlots:
        PlotHelper.savePlot(figEtaINScan, "etaINScan.{}".format(savePlots))
    #}}}

    #{{{etaE plot
    for nr, key in enumerate(sortedEtaEKeys):
        if  B0 == B0s[0]:
            label = "${}$".format(key)
        else:
            label = None
        axEtaENScan.plot(n0s,\
                         etaE[key],\
                         color = colors[nr],\
                         ls    = ls,\
                         label = label)

    # Set axis
    axEtaENScan.set_yscale("log")
    setNScanAxes(figEtaENScan, axEtaENScan, r"$\eta_e/m_i n \rho_s c_s$", Te0EV)

    if savePlots:
        PlotHelper.savePlot(figEtaENScan, "etaENScan.{}".format(savePlots))
    #}}}
    #}}}

    #{{{Te scan
    # Fixed n0
    n0 = 1e17
    # Convert to J
    Te0s = Te0sEV*cst.eV

    # Calculate normalization variables
    cS   = p.calcCS(Te0, Ti0, N)
    rhoS = p.calcRhoS(cS, omCI)

    # Caclulate collisionalities
    coulombLog = p.calcCoulombLog(Te0s, n0)
    vThE = p.calcVThE(Te0s)
    vThI = p.calcVThI(Ti0)
    nuEI = p.calcNuEI(n0, coulombLog, vThE)
    nuII = p.calcNuII(n0, coulombLog, vThI)
    tauI = 1.0/(sqrt(2.0)*nuII)
    tauE = 1.0/nuEI

    # Caclulate ion viscosities
    # Make an array in order to get correct dimensions
    omCIArray = np.zeros(tauI.shape)
    omCIArray.fill(omCI)
    etaI = p.calcEtaI(n0, Ti0, tauI, omCIArray)
    sortedEtaIKeys = list(etaI.keys())
    sortedEtaIKeys.sort()

    # Caclulate electron viscosities
    etaE = p.calcEtaE(n0s, Te0s, tauE, omCE)
    sortedEtaEKeys = list(etaE.keys())
    sortedEtaEKeys.sort()

    # Normalization
    nuEI /= omCI
    nuII /= omCI

    for key in sortedEtaIKeys:
        etaI[key] /= (mi*n0*rhoS*cS)

    for key in sortedEtaEKeys:
        etaE[key] /= (mi*n0*rhoS*cS)

    # Normalization
    nuEI /= omCI
    nuII /= omCI

    #{{{nuEI plot
    axNuEITScan.plot(Te0sEV,\
                     nuEI,\
                     color=colors[0],\
                     ls   = ls,\
                     label="$B_0={}$".format(B0))

    # Set axis
    setTScanAxes(figNuEITScan, axNuEITScan, r"$\nu_{ei}/\omega_{ci}$", n0)

    if savePlots:
        PlotHelper.savePlot(figNuEITScan, "nuEITScan.{}".format(savePlots))
    #}}}

    #{{{etaI plot
    for nr, key in enumerate(sortedEtaIKeys):
        if  B0 == B0s[0]:
            label = "${}$".format(key)
        else:
            label = None
        axEtaITScan.plot(Te0sEV,\
                         etaI[key],\
                         color = colors[nr],\
                         ls    = ls,\
                         label = label)

    # Set axis
    axEtaITScan.set_yscale("log")
    setTScanAxes(figEtaITScan, axEtaITScan, r"$\eta_i/m_i n \rho_s c_s$", n0)

    if savePlots:
        PlotHelper.savePlot(figEtaITScan, "etaITScan.{}".format(savePlots))
    #}}}

    #{{{etaE plot
    for nr, key in enumerate(sortedEtaEKeys):
        if  B0 == B0s[0]:
            label = "${}$".format(key)
        else:
            label = None
        axEtaETScan.plot(Te0sEV,\
                         etaE[key],\
                         color = colors[nr],\
                         ls    = ls,\
                         label = label)

    # Set axis
    axEtaETScan.set_yscale("log")
    setTScanAxes(figEtaETScan, axEtaETScan, r"$\eta_e/m_i n \rho_s c_s$", n0)

    if savePlots:
        PlotHelper.savePlot(figEtaETScan, "etaETScan.{}".format(savePlots))
    #}}}
    #}}}

if showPlots:
    # Save the plot
    plt.show()

plt.close("all")
