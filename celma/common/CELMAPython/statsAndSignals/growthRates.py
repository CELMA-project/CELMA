#!/usr/bin/env python

"""
Contains functions to calculate and plot the angular frequency and
growth rate of a time trace of a spatial FFT
"""

from ..plotHelpers import PlotHelper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec

#{{{linRegOfExp
def linRegOfExp(x,y):
    #{{{docstring
    """
    Calculates the gradient of an exponential function using linear
    regression.

    NOTE: There seem to be a lot of confusion of how to estimate
          uncertainties. To be consise, we will here follow
          "An introduction to error analysis" by Tylor, J.R, and
          not use what is used in polyfit or stats.linregress.

    NOTE: We are fitting an exponential function. The uncertainties in the
          gradient is strictly speaking NOT given in (8.17).
          We should use the weigthed fit, using (8.34) and the
          equation of p 199 to  estimate the uncertainties in
          sigma_b.
          However, we are here only interested to get the spread in the
          data of the straigth line (the standard deviation of the
          straigth line), and equation (8.17) will be used. If we were
          to calculate the spread in the exponential data (assuming that
          the uncertainties were constant), we would have

            >>> # Calculation of weigths
            >>> # We have that
            >>> # lnY = ln(y).
            >>> # We have assumed that the uncertainties in y is equally
            >>> # uncertain. From error propagation, we then get.
            >>> sigmaY = sigmaLnY/y
            >>> w      = 1/sigmaY**2
            >>> # Regression of exponential
            >>> DeltaExp = sum(w)*sum(w*x**2) - (sum(w*x))**2
            >>> BExp     = (sum(w)*sum(w*lnY*x) - sum(w*lnY)*sum(w*x))/DeltaExp
            >>> sigmaB   = np.sqrt(sum(w)/DeltaExp)

          Note also that there is not a huge difference between B and BExp

    Parameters
    ----------
    x : array
        The x data. Should be without measuring uncertainties.
    y : array
        The y data. It will be assumed that the uncertainties in y is
        equally uncertain. This means that the uncertainties in ln(y)
        will scale as sigmaLnY/y

    Returns
    -------
    B : float
        The exponential growth rate
    sigmaB : float
        The uncertainties in BExp (see notes above)
    """
    #}}}

    # Take the logarithm
    lnY = np.log(y)
    # Calculation of sigmaLnY
    N        = len(x)
    Delta    = N*sum(x**2) - (sum(x))**2
    A        = (sum(x**2)*sum(lnY) - sum(x)*sum(x*lnY))/Delta
    B        = (N*sum(lnY*x) - sum(x)*sum(lnY))/Delta
    sigmaLnY = np.sqrt((1/(N-2))*sum((lnY - A - B*x)**2))
    sigmaB   = sigmaLnY*np.sqrt(N/Delta)

    return B, sigmaB
#}}}

#{{{calcGrowthRate
def calcGrowthRate(modes, time, maxMode = 7, diagnose=False):
    #{{{docstring
    r"""
    Calculates the angular frequency and growth rate of a time trace of
    a spatial FFT.

    Analytic results are often obtained from linear theory, and are
    often assuming that the perturbation acts like plane waves, that is

    \tilde{f}_1 = A\cdot\exp(-i[\mathbf{k}\cdot\mathbf{x} + \omega t]) + cc.

    In that case, the growth rate would be interpreted as

    \Im{\omega}

    and the angular frequency would be interpreted as

    \Re{\omega}

    To compare with linear theory, the timetrace before the onset of
    turbulence must be considered, where the perturbations are small so
    that mode coupling can be negleted and the system is basically
    acting linearly.

    This routine assumes exponential growth rates (i.e. the logarithm
    will be taken of a given mode, and the growth rate is found from the
    linear regression of the line).

    In order to find straight lines, the signal is first binned, and the
    gradient of each bin is found using linear regression. If the
    difference in the gradient is less than 5% for successive bins, the
    bins will be registered as a hit. A straigth line in the logarithmic
    signal is defined to be where there has been at least 4 hits.

    The angular velocity is found from the phase difference of subsequent
    times within the time the signal is registered as a straight line.

    Parameters
    ----------
    modes : 2d-array
        The time traces of the FFTed signal. The first index represents
        the time index and the second index represents the mode index
    time : array
        Time corresponding to the modes
    maxMode : int
        How many modes to investigate
    diagnose : bool
        Will print extra diagnostics if True

    Returns
    -------
    results : [dict | None]
        A dictionary of dictionary. The top level key is the mode number
        and the bottom level keys are
            * angFreq       - The angular frequency
            * angFreqStd    - The standard deviation of the angular frequency
            * growthRate    - The growth rate
            * growthRateStd - The standard deviation of the growth rate
            * startTime     - Start time of the linear segment
            * endTime       - End time of the linear segment
        None is returned if exception occured
    """
    #}}}

    #{{{ NOTE: We are dealing with a real signal:
    #          As the fourier transform breaks the signal up in cisoids
    #          there will be one part of the signal in the positive
    #          rotating ciscoid and one in the negative (negative
    #          frequencies) for a given mode number. We need to take
    #          into account both in order to calculate the amplitude. As
    #          the signal is real only one of the phase sifts are
    #          needed. Notice that for a real signal the imaginary part
    #          occurs as a complex conjugate pair
    # http://dsp.stackexchange.com/questions/431/what-is-the-physical-significance-of-negative-frequencies?noredirect=1&lq=1
    # http://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
    #}}}

    # Number of points
    N = modes.shape[1]

    # This is part of the algorithm used to detect the straigth line in
    # the logarithm of the signal
    # Bin the time axis
    binSize = 20
    bins = np.arange(0,len(modes[:,0]),binSize)

    if bins.shape[0] == 0:
        message = ("{0}{1}WARNING: "\
                   "No proper bins could be made, returning None{1}{0}")
        print(message.format("\n"*2, "!"*5))
        return None

    # Finding mode number
    results = {}

    # Start on 1 as we have neglegted the DC mode
    # +1 as we skip the DC mode
    # +1 as range excludes the last point
    for mNr in range(1, maxMode+2):
        # Place holders for the growth rates and the mean square error
        growthRates  = []
        startIndices = []

        # Loop over the bins, start from 1 as we will index b-1
        for bNr, b in enumerate(range(1, len(bins))):
            # Find the growth rate of the current bin
            startIndex = bins[b-1]
            endIndex   = bins[b]
            curTime    = time[startIndex: endIndex]
            # Magnitude of the signal
            # https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition
            modeMag = (np.abs(modes[startIndex: endIndex,  mNr]) +\
                       np.abs(modes[startIndex: endIndex, -mNr]))/N

            growthRate, _ = linRegOfExp(curTime, modeMag)
            growthRates .append(growthRate)
            startIndices.append(startIndex)

            if diagnose:
                # Print diagnostics
                message = ("mode = {:<8} startIndex={:<8} endIndex={:<8} "
                           "growthRate={:<+8.3f} ")
                print(message.format(mNr, startIndex, endIndex, growthRate))

        #{{{Find growth rates from definition of straight segments in the plot
        # Initialize the previous growth rate
        prevRate = growthRates[0]

        # The growth rates and the corresponding standard deviation is found by
        # the criteria that the growth rate should stay the same ober at least
        # 4 bins. If the criterion is met, a new linear regression will be
        # taken from the max and min of the time indices.

        # Place holders (used in case there are several hits)
        finalGrowthRates  = []
        sigmaBs           = []
        finalStartIndices = []
        finalEndIndices   = []
        startIndex = 0
        hits       = 0

        # Looping from 1 as we already have the previous growth rate
        for gNr, growthRate in enumerate(growthRates[1:]):
            # If the growth rate is within 5%
            if growthRate*0.95 <= prevRate <= growthRate*1.05:
                if hits == 0:
                    # Include the index of the previous (not subtracting
                    # with 1 as gNr lags the growthRates by 1)
                    startIndex = startIndices[gNr]
                    # Update hits counter
                    hits += 1
                # Update hits counter
                hits += 1
            else:
                if hits >= 4:
                    endIndex = bins[int(startIndex/binSize) + hits]
                    curTime  = time[startIndex: endIndex]
                    modeMag  = \
                        (np.abs(modes[startIndex: endIndex,  mNr]) +\
                         np.abs(modes[startIndex: endIndex, -mNr]))/N
                    growthRate, sigmaB = linRegOfExp(curTime, modeMag)
                    finalGrowthRates.append(growthRate)
                    sigmaBs    .append(sigmaB)
                    # We will currently use indices for the start times
                    # as these are easier to deal with when finding the
                    # angular velocity. The indices will be converted to
                    # actual time in the end
                    finalStartIndices.append(startIndex)
                    finalEndIndices  .append(endIndex)
                # Reset the counter
                hits = 0
            prevRate = growthRate

        # If the last element in the growthRates gave a hit
        if hits >= 4:
            endIndex = bins[int(startIndex/binSize) + hits]
            curTime  = time[startIndex: endIndex]
            modeMag  = \
                (np.abs(modes[startIndex: endIndex,  mNr]) +\
                 np.abs(modes[startIndex: endIndex, -mNr]))/N
            growthRate, sigmaB = linRegOfExp(curTime, modeMag)
            finalGrowthRates .append(growthRate)
            sigmaBs          .append(sigmaB)
            finalStartIndices.append(startIndex)
            finalEndIndices  .append(endIndex)

        if len(finalGrowthRates) > 0:
            # Select the first growth rate
            finalStartIndex = finalStartIndices[0]
            finalEndIndex   = finalEndIndices  [0]
            results[mNr]    = {\
                    "growthRate"    : finalGrowthRates[0              ],\
                    "growthRateStd" : sigmaBs         [0              ],\
                    "startTime"     : time            [finalStartIndex],\
                    "endTime"       : time            [finalEndIndex  ],\
            }

            if len(finalGrowthRates) > 1:
                message = ("{0}{1}WARNING: "\
                           "Found {2} different linear segments. "\
                           "Selecting the first in the key 'growthRates'. "\
                           "Storing all growth rates in the key "\
                           "'allGrowthRates'{1}{0}")
                results[mNr]["allGrowthRates"] = finalGrowthRates
                print(message.format("\n"*2, "!"*5, len(finalGrowthRates)))
        else:
            results[mNr] = None
        #}}}

        #{{{Finding the real part
        # If no growth rate was found
        if results[mNr] is None:
            continue
        deltaT = time[1] - time[0]
        # Create the place holder for the angular frequency
        angularFreq = np.zeros(finalEndIndex-finalStartIndex)
        # Loop over the indices between start and stop
        # finalStartIndex+1 as we will index i-1
        # finalEndIndex+1 to include the last point
        for nr, i in enumerate(range(finalStartIndex+1, finalEndIndex+1)):
            # The phase shift is found from atan2
            # http://dsp.stackexchange.com/questions/23994/meaning-of-real-and-imaginary-part-of-fourier-transform-of-a-signal
            # atan2 in [-pi, pi]
            prevPhaseShift =\
                np.arctan2(modes[i-1, mNr].imag, modes[i-1, mNr].real)
            curPhaseShift  =\
                np.arctan2(modes[i  , mNr].imag, modes[i  , mNr].real)
            # phaseShiftDiff in [0, 2*pi]
            phaseShiftDiff = prevPhaseShift - curPhaseShift

            if curPhaseShift*prevPhaseShift < 0\
               and abs(curPhaseShift) + abs(prevPhaseShift) > np.pi:
                if curPhaseShift < 0:
                    if diagnose:
                        print("Phase shift shifted from pi to -pi")
                    # We are going from pi to -pi
                    # In order to avoid the discontinuity, we turn
                    # curPhaseShift and prevPhaseShift to the opposite
                    # quadrants
                    tempCurPhase   = np.pi + curPhaseShift
                    tempPrevPhase  = np.pi - prevPhaseShift
                    phaseShiftDiff = -(tempCurPhase + tempPrevPhase)
                else:
                    if diagnose:
                        print("Phase shift shifted from -pi to pi")
                    # We are going from -pi to pi
                    # In order to avoid the discontinuity, we turn
                    #curPhaseShift and prevPhaseShift to the opposite quadrants
                    tempCurPhase   = np.pi - curPhaseShift
                    tempPrevPhase  = np.pi + prevPhaseShift
                    phaseShiftDiff = tempCurPhase + tempPrevPhase

            # The angular speed (angular frequency) has units rad/s.
            # Remember that if angularFreq*t = 2*pi the perturbation has
            # revolved one time
            angularFreq[nr] = phaseShiftDiff/deltaT

            if diagnose:
                # Print diagnostics
                message = ("mode = {:<8} index={:<8} angFreq={:<+8.3f} "\
                           "phase={:<+8.3f} phaseDiff={:<+8.3f}")
                print(message.format(mNr, i, angularFreq[nr],\
                                     curPhaseShift, phaseShiftDiff))

        # Calculate the mean and the spread (the standard deviation), note
        # that numpys std misses a minus 1 in the denominator, but as N is
        # high, this is negligible
        # NOTE: A negative angular frequency means that the
        #       perturbations are moving in the negative theta
        #       direction.
        results[mNr]["angFreq"]    = angularFreq.mean()
        results[mNr]["angFreqStd"] = angularFreq.std()
        #}}}

    return results
#}}}

#{{{plotGrowthRates
def plotGrowthRates(df):
    #{{{docstring
    """
    Plots the growth rates

    Parameters
    ----------
    df : DataFrame
        The data frame to be plotted. The data frame must be orgnized
        such that:
        1. The x axis will be given by splitting df.columns by "=" and take
           the value right og the "="
        2. The index only has two levels
        3. The plot name is given at the first level of df.index
        4. The growth rate, angular frequency and their standar
           deviation is given in the innermost level of df
    """
    #}}}

    # FIXME: T ALWAYS NORMALIZED WITH OMCI
    # FIXME: Other normalizations
    #{{{ FIXME: Consider make this a class
    # Dict used as an equivalent to C++ maps
    mapToPltText = {"modeNr":"$Mode number$",\
                    "B0"    :"$B_0$"        ,\
                    "Te0"   :"$T_e$"        ,\
                    "nn"    :"$n_n$"        ,\
                   }

    errorbarOptions = {"color"     :"k",\
                       "fmt"       :"o",\
                       "fmt"       :"o",\
                       "markersize":10 ,\
                       "ecolor"    :"k",\
                       "capsize"   :7  ,\
                       "elinewidth":3  ,\
                       }

    figSize = (9,12)
    markeredgewidth = 3

    showPlot = False
    savePlot = True
    extension = "png"

    All = slice(None)
    #}}}


    # Loop through the figures
    for plotLabel in df.index.levels[0]:
        fig = plt.figure(figsize = figSize)
        gs  = GridSpec(nrows=2, ncols=1)

        imAx   = fig.add_subplot(gs[0])
        realAx = fig.add_subplot(gs[1], sharex=imAx)

        # We can now access unsing the loc method
        # The syntax is
        # df.loc[("indexLevel1", "indexLevel2", ...),\
        #        ("columnsLevel1", "columnsLevel2", ...)]
        # http://pandas.pydata.org/pandas-docs/stable/cookbook.html#cookbook-selection
        xAxis   = [float(txt.split("=")[1])\
                   for txt in df.loc[(plotLabel, "growthRate"), All].index]
        yAxisIm = df.loc[(plotLabel, "growthRate"), All].values
        yErrIm  = df.loc[(plotLabel, "growthRateStd"), All].values
        yAxisRe = df.loc[(plotLabel, "angularFreq"), All].values
        yErrRe  = df.loc[(plotLabel, "angularFreq"), All].values

# FIXME: Instead of errorbar use a partial function?
        (_, caps, _) = imAx.errorbar(xAxis,\
                                     yAxisIm,\
                                     yerr=yErrIm,\
                                     **errorbarOptions)

        # Set capsize
        for cap in caps:
            cap.set_markeredgewidth(markeredgewidth)

        # Angular frequencies
        (_, caps, _) = realAx.errorbar(xAxis,\
                                       yAxisRe,\
                                       yerr=yErrRe,\
                                       **errorbarOptions)

        # Set capsize
        for cap in caps:
            cap.set_markeredgewidth(markeredgewidth)

        # Add 10% margins for readability
        imAx  .margins(x=0.1, y=0.1)
        realAx.margins(x=0.1, y=0.1)

        # FIXME: Rotation
        PlotHelper.makePlotPretty(imAx,   yprune = "both", rotation = None)
        PlotHelper.makePlotPretty(realAx, yprune = "both", rotation = None)

        fig.suptitle(mapToPltText[plotLabel.split("=")[0]])
        imAx  .set_ylabel("$\omega_I$ $[]$")
        realAx.set_ylabel("$\omega_R$ $[]$")

        imAx.tick_params(labelbottom="off")
        xlabel = mapToPltText[\
                    df.loc[(plotLabel, "growthRate"), All].index[0].\
                    split("=")[0]\
                 ]
        realAx.set_xlabel(xlabel)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0)

        # Sort the xData and use it as ticks
        xAxis.sort()
        realAx.xaxis.set_ticks(xAxis)

        if savePlot:
            PlotHelper.savePlot(fig,\
                "growthrate_{}_{}.{}".\
                format(plotLabel, xlabel, extension))

        if showPlot:
            plt.show()
#}}}
