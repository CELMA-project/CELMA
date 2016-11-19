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
from matplotlib.ticker import MaxNLocator

showPlots = False
savePlots = ".pdf"

#{{{Common for simple plots
# Set the plotting style
title_size = 30
title_height = 1.02
plt.rc("font", size = 30)
plt.rc("axes", labelsize = 25, titlesize = title_size)
plt.rc("xtick", labelsize = 25)
plt.rc("ytick", labelsize = 25)
plt.rc("legend", fontsize = 30)
plt.rc("lines", linewidth = 3.0)
plt.rc("lines", markersize = 20.0)
plt_size = (10, 7)
fig_no = 1
# Try to make a figure with the current backend
try:
    fig = plt.figure(fig_no, figsize = plt_size)
except:
    # Switch if a backend needs the display
    plt.switch_backend('Agg')
    fig = plt.figure(fig_no, figsize = plt_size)
#}}}}

#{{{setNScanAxes
def setNScanAxes(fig, ax, yLabel, Te0EV):
    ax.set_xlabel(r"$n [\mathrm{m}^{-3}]$")
    ax.set_ylabel(yLabel)
    title = r'$Te=${} $[\mathrm{{eV}}]$'.format(plotNumberFormatter(Te0EV))
    fig.suptitle(title, y=title_height)
    makeAxPretty(ax)

    # Includes the xlabel if outside
    fig.tight_layout(rect=[0.001,0.001,1,1])
#}}}

#{{{setTScanAxes
def setTScanAxes(fig, ax, yLabel, n0):
    # Set axis label
    ax.set_xlabel(r"$Te [\mathrm{eV}]$")
    ax.set_ylabel(yLabel)
    title = r'$n=${} $[\mathrm{{m}}^{{-3}}]$'.format(plotNumberFormatter(n0))
    fig.suptitle(title, y=title_height)
    makeAxPretty(ax)

    # Includes the xlabel if outside
    fig.tight_layout(rect=[0.001,0.001,1,1])
#}}}

#{{{plotNumberFormatter
def plotNumberFormatter(val):
    """
    Formatting numbers in the plot

    Input
    val - The value
    """

    tickString = '${:.3g}'.format(val)
    if "e+" in tickString:
        tickString = tickString.replace('e+', r'\cdot 10^{')
        tickString += '}$'
    elif "e-" in tickString:
        tickString = tickString.replace('e-', r'\cdot 10^{-')
        tickString += '}$'
    else:
        tickString += '$'

    return tickString
#}}}

#{{{makeAxPretty
def makeAxPretty(ax):
    # Reset the tick labels
    tickLabels = [plotNumberFormatter(el) for el in ax.get_xticks().tolist()]
    ax.set_xticklabels(tickLabels, rotation='vertical')
    tickLabels = [plotNumberFormatter(el) for el in ax.get_yticks().tolist()]
    ax.set_yticklabels(tickLabels)
    # Plot the legend
    leg = ax.legend(loc="best", fancybox = True, numpoints=1)
    leg.get_frame().set_alpha(0.5)
    # Plot the grid
    ax.grid()
    # Make sure no collision between the ticks
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
#}}}

#{{{savePlot
def savePlot(fig, name):
    fig.savefig(name + savePlots,\
                transparent = True    ,\
                bbox_inches = 'tight' ,\
                pad_inches  = 0       ,\
                )
#}}}

mi = cst.m_p
p = Parameters(cst.m_p)

B0s   = [1, 1e-1, 1e-2]
Ti0EV = 7.0
n0s    = np.linspace(1e16, 1e18, 200)
Te0sEV = np.linspace(7, 21, 200)

B0LineType = ['solid', 'dashed', 'dotted', 'dash_dot']
# One color for each eta
colors = cm.rainbow(np.linspace(0, 1, 5))

# Make axes
figNuEINScan = plt.figure(0, figsize = plt_size)
figNuEITScan = plt.figure(1, figsize = plt_size)
figEtaINScan = plt.figure(2, figsize = plt_size)
figEtaITScan = plt.figure(3, figsize = plt_size)
figEtaENScan = plt.figure(4, figsize = plt_size)
figEtaETScan = plt.figure(5, figsize = plt_size)
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
        savePlot(figNuEINScan, "niEINScan")
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
    axEtaINScan.set_yscale('log')
    setNScanAxes(figEtaINScan, axEtaINScan, r"$\eta_i/m_i n \rho_S c_S$", Te0EV)

    if savePlots:
        savePlot(figEtaINScan, "etaINScan")
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
    axEtaENScan.set_yscale('log')
    setNScanAxes(figEtaENScan, axEtaENScan, r"$\eta_e/m_i n \rho_S c_S$", Te0EV)

    if savePlots:
        savePlot(figEtaENScan, "etaENScan")
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
    etaI = p.calcEtaI(n0, Ti0, tauI, omCI)
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
        savePlot(figNuEITScan, "nuEITScan")
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
    axEtaITScan.set_yscale('log')
    setTScanAxes(figEtaITScan, axEtaITScan, r"$\eta_i/m_i n \rho_S c_S$", n0)

    if savePlots:
        savePlot(figEtaITScan, "etaITScan")
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
    axEtaETScan.set_yscale('log')
    setTScanAxes(figEtaETScan, axEtaETScan, r"$\eta_e/m_i n \rho_S c_S$", n0)

    if savePlots:
        savePlot(figEtaETScan, "etaETScan")
    #}}}
    #}}}

if showPlots:
    # Save the plot
    plt.show()

plt.close("all")
