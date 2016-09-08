#!/usr/bin/env python

import numpy as np
from numpy import log, sqrt
import scipy.constants as cst
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

#{{{calcCoulombLog
def calcCoulombLog(Te0, n0):
    coulombLog = log(12.0*pi*(e0*Te0)**(3/2)/(sqrt(n0)*q**3))
    return coulombLog
#}}}

#{{{calcVThE
def calcVThE(Te0):
    vThE = sqrt(2.0*Te0/me)
    return vThE
#}}}

#{{{calcVThI
def calcVThI(Ti0):
    vThE = sqrt(2.0*Ti0/mi)
    return vThE
#}}}

#{{{calcNuEI
def calcNuEI(ne0, coulombLog, vThE):
    nuEI=1.0/(4.0*pi)*(q**4.0*ne0/(e0**2.0*me**2.0))*coulombLog*(1.0/vThE**3.0)
    return nuEI
#}}}

#{{{calcNuII
def calcNuII(ne0, coulombLog, vThI):
    nuII=1.0/(2.0*pi)*(q**4.0*ne0/(e0**2.0*mi**2.0))*coulombLog*(1.0/vThI**3.0)
    return nuII
#}}}

#{{{calcEtaI
def calcEtaI(n, Ti, tauI, omCI):
    etaI = {}

    etaI['\eta_{{i,0}}'] = 0.96*n*Ti*tauI
    etaI['\eta_{{i,1}}'] = 3.0*n*Ti/(10.0*(omCI**2.0)*tauI)
    etaI['\eta_{{i,2}}'] = 4.0*etaI['\eta_{{i,1}}']
    etaI['\eta_{{i,3}}'] = n*Ti/(2.0*omCI)
    etaI['\eta_{{i,4}}'] = 2.0*etaI['\eta_{{i,3}}']

    return etaI
#}}}

#{{{calcEtaE
def calcEtaE(n, Te, tauE, omCe):
    etaE = {}

    etaE['\eta_{{e,0}}'] = 0.73*n*Te*tauE
    etaE['\eta_{{e,1}}'] = 0.51*n*Te/(omCE**2.0*tauE)
    etaE['\eta_{{e,2}}'] = 4.0*etaE['\eta_{{e,1}}']
    etaE['\eta_{{e,3}}'] = n*Te/(2.0*omCE)
    etaE['\eta_{{e,4}}'] = 2.0*etaE['\eta_{{e,3}}']

    return etaE
#}}}

# Aliases
q   = cst.e
e0  = cst.epsilon_0
me  = cst.m_e
mi  = cst.m_p
pi  = np.pi

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
    omCI = q*B0/mi
    omCE = q*B0/me
    N    = 3.0

    # Convert to J
    Ti0 = Ti0EV*cst.eV

    #{{{n scan
    # Fixed Te0
    Te0EV = 14.0
    # Convert to J
    Te0 = Te0EV*cst.eV

    # Calculate normalization variables
    cS   = sqrt((Te0+((N+2.0)/N)*Ti0)/mi)
    rhoS = cS/omCI

    # Caclulate collisionalities
    coulombLog = calcCoulombLog(Te0, n0s)
    vThE = calcVThE(Te0)
    vThI = calcVThI(Ti0)
    nuEI = calcNuEI(n0s, coulombLog, vThE)
    nuII = calcNuII(n0s, coulombLog, vThI)
    tauI = 1.0/(sqrt(2.0)*nuII)
    tauE = 1.0/nuEI

    # Caclulate ion viscosities
    etaI = calcEtaI(n0s, Ti0, tauI, omCI)
    sortedEtaIKeys = list(etaI.keys())
    sortedEtaIKeys.sort()

    # Caclulate electron viscosities
    etaE = calcEtaE(n0s, Te0, tauE, omCE)
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
    cS   = sqrt((Te0s+((N+2.0)/N)*Ti0)/mi)
    rhoS = cS/omCI

    # Caclulate collisionalities
    coulombLog = calcCoulombLog(Te0s, n0)
    vThE = calcVThE(Te0s)
    vThI = calcVThI(Ti0)
    nuEI = calcNuEI(n0, coulombLog, vThE)
    nuII = calcNuII(n0, coulombLog, vThI)
    tauI = 1.0/(sqrt(2.0)*nuII)
    tauE = 1.0/nuEI

    # Caclulate ion viscosities
    etaI = calcEtaI(n0, Ti0, tauI, omCI)
    sortedEtaIKeys = list(etaI.keys())
    sortedEtaIKeys.sort()

    # Caclulate electron viscosities
    etaE = calcEtaE(n0s, Te0s, tauE, omCE)
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
