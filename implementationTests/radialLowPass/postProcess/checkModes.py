#!/usr/bin/env python

"""Check the modes of the run"""

from functools import partial
from boutdata import collect
import numpy as np
import matplotlib.pylab as plt

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../../celma/common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotHelpers import PlotHelper

# Set the plot style for all plots
titleSize = 30
plt.rc("font",   size      = 30)
plt.rc("axes",   labelsize = 25, titlesize = titleSize)
plt.rc("xtick",  labelsize = 25)
plt.rc("ytick",  labelsize = 25)
plt.rc("legend", fontsize  = 20)
plt.rc("lines",  linewidth = 2)


def checkModes(path):
    defaultCollect = partial(collect,\
                        path=path, xguards=False, yguards=False, info=False)
    unfiltered = defaultCollect("unfiltered")
    filtered   = defaultCollect("filtered"  )
    Lx         = defaultCollect("Lx"        )
    dx         = defaultCollect("dx"        )

    maxRhoInd = unfiltered.shape[0]
    # NOTE: maxRhoInd counts from 0
    indices = (0, int(maxRhoInd/2)-1, maxRhoInd-1)

    filteredDict   = {}
    unfilteredDict = {}
    fftFiltered    = {}
    fftUnfiltered  = {}

    # Print stats
    kMax = np.floor((filtered.shape[-1]/2)*(2/3))
    print(("\n\nMax allowed mode number on outer circumference = Nyquist mode "
            "* Orszag = floor((nz/2)*(2/3)) = floor(({}/2)*(2/3)) = {} ")\
            .format(filtered.shape[-1], kMax)
         )
    lambdaMin = 2*np.pi*(Lx-dx[0,0]/2)/kMax
    print(("\nThis means that the minimum wavelength"
           " lambdaMin = outer circmference/kMax ="
           " 2*pi*(rho-dx/2)/kMax = 2*pi*({}-{}/2)/{} = {}")\
           .format(Lx, dx[0,0], kMax, lambdaMin)
         )
    print("\nThis means that kCurX = floor(cur circumference/lambdaMin)")
    for ind in range(filtered.shape[0]):
        # The first inner point is 0.5*dx from the origin
        x       = dx[0,0]*(ind + 0.5)
        curCirc = 2*np.pi*x
        kCurXClean = curCirc/lambdaMin
        kCurX      = np.floor(curCirc/lambdaMin)
        # Find first index where the value is close to 0
        firstZero\
            = np.where(\
                np.isclose(\
                    np.fft.fft(\
                        filtered[ind, 0, :].flatten().real), 0.0))[0][0]

        print(("Ind={:<3} => x={:<7.2f}=> cur circumference={:<7.2f}"
               " => kCurX={:<2d}, floored from {:<9.6f}. Last value found"
               " on mode {}")\
               .format(ind, x, curCirc, int(kCurX), kCurXClean, firstZero-1)
             )

    for ind in indices:
        unfilteredDict[ind] = unfiltered[ind, 0, :].flatten()
        filteredDict  [ind] = filtered[ind, 0, :]  .flatten()
        fftUnfiltered [ind] = np.fft.fft(unfilteredDict[ind])
        fftFiltered   [ind] = np.fft.fft(filteredDict  [ind])

    fig = plt.figure()
    axUnfilter    = fig.add_subplot(221)
    axFilter      = fig.add_subplot(223)
    axFFTUnfilter = fig.add_subplot(222)
    axFFTFilter   = fig.add_subplot(224)

    axUnfilter   .set_xlabel("Unfiltered")
    axFilter     .set_xlabel("Filtered")
    axFFTUnfilter.set_xlabel("FFT Unfiltered")
    axFFTFilter  .set_xlabel("FFT Filtered")

    for ind in indices:
        axUnfilter   .plot(unfilteredDict[ind].real,"o",alpha=0.5,label="Index {}".format(ind))
        axFilter     .plot(filteredDict  [ind].real,"o",alpha=0.5,label="Index {}".format(ind))
        axFFTUnfilter.plot(fftUnfiltered [ind].real,"o",alpha=0.5,label="Index {}".format(ind))
        axFFTFilter  .plot(fftFiltered   [ind].real,"o",alpha=0.5,label="Index {}".format(ind))

    PlotHelper.makePlotPretty(axUnfilter)
    PlotHelper.makePlotPretty(axFilter)
    PlotHelper.makePlotPretty(axFFTUnfilter)
    PlotHelper.makePlotPretty(axFFTFilter)

    plt.show()
