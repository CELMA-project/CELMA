#!/usr/bin/env python

"""
Creation of the parallel current balance.
"""

import pickle
import matplotlib.pylab as plt
import numpy as np
import scipy.constants as cst
from boututils.options import BOUTOptions

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../../../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.plotHelpers import (SizeMaker,\
                                 PlotHelper,\
                                 plotNumberFormatter,\
                                 seqCMap3)

# Collect from fieldScan
path = "../../CSDXMagFieldScanAr/visualizationNormalized/B0_0.06/field1D/jPar-parallel-1D-0.pickle"
with open(path, "rb") as f:
    fig = pickle.load(f)

axes = fig.get_axes()

for ax in axes[:-1]:
    fig.delaxes(ax)

oldAx = axes[-1]

# Convert to physical
inputFileOpts = BOUTOptions("../../CSDXNyScan")
n0  = eval(inputFileOpts.input["n0"])
Te0 = eval(inputFileOpts.input["Te0"])*cst.e
B0  = eval(inputFileOpts.input["B0"])
if inputFileOpts.input["gas"] != "Ar":
    raise ValueError("This routine is only made for Ar gas")

mi  = 39.948*cst.u
e   = cst.e

omCI = e*B0/mi
cS   = np.sqrt(Te0/mi)
rhoS = cS/omCI

# From normalization of parallel current equation
factor = omCI*n0*e*cS

# Obtain the suptitle
suptitle = fig.texts[0].get_text()

# Make new ax to plot to
fig, newAx = plt.subplots(figsize = SizeMaker.standard(w=4, a=0.5))

maxmin = []
for line in oldAx.get_lines():
    x, y = line.get_data()
    x *= rhoS
    y *= factor
    ny = len(x)
    maxmin.append((np.max(y), np.min(y)))
    newAx.plot(x,y, color=line.get_color())

# Set the texts
maxInd = np.argmax(tuple(curVal[0] for curVal in maxmin))
minInd = np.argmin(tuple(curVal[1] for curVal in maxmin))

boltzmann   = oldAx.get_lines()[minInd]
resistivity = oldAx.get_lines()[maxInd]

boltzmannTxt = r"$\mu n\partial_{\parallel}\left(T_e \ln(n) - \phi \right)$"
newAx.text(0.8, -1.8e8, boltzmannTxt, color=boltzmann.get_color(),\
           va="center", ha="center", size="x-large",\
           bbox={"facecolor":"white",\
                 "edgecolor":"none" ,\
                 "boxstyle" :"round",\
                 "alpha"    :0.9},\
           )

resTxt = r"$-0.51\nu_{ei}j_\parallel$"
newAx.text(0.7, 1.8e8, resTxt, color=resistivity.get_color(),\
           va="center", ha="center", size="x-large",\
           bbox={"facecolor":"white",\
                 "edgecolor":"none" ,\
                 "boxstyle" :"round",\
                 "alpha"    :0.9},\
           )

ddtTxt = r"$\partial_t j_\parallel$"
newAx.text(1.4, 0, ddtTxt, color="k",\
           va="center", ha="center", size="x-large",\
           bbox={"facecolor":"white",\
                 "edgecolor":"none" ,\
                 "boxstyle" :"round",\
                 "alpha"    :0.9},\
           )

newAx.set_xlabel(r"$z\;[\mathrm{m}]$")
newAx.set_ylabel(r"$[\mathrm{C}\mathrm{m}^{-2}\mathrm{s}^{-2}]$")

# Create new title
theSplit = suptitle.split("=")
rho = float(theSplit[1].split("$")[2])*rhoS
t   = eval(theSplit[-1].split("$")[2].\
           replace("{","").replace("}","").replace("\\cdot 10^","e"))/omCI
newAx.set_title((r"$\rho={}\;\mathrm{{m}}\quad"
                 r"\theta=0^{{\circ}}\quad t={}\;\mathrm{{s}}$").\
        format(plotNumberFormatter(rho,0).replace("$",""),\
               plotNumberFormatter(t,0)  .replace("$","")\
        ))

PlotHelper.makePlotPretty(newAx, legend=False)
PlotHelper.savePlot(fig, "jParBalanceNy{}.pdf".format(ny))
