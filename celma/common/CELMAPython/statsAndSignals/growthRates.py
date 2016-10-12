#!/usr/bin/env python

"""
Contains functions to calculate and plot the angular frequency and
growth rate of a time trace of a spatial FFT
"""

import numpy as np

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
def calcGrowthRate(modes, time, maxMode = 7):
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
    difference in the gradient is less than 1% for successive bins, the
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

    Returns
    -------
    results : dict
        A dictionary of dictionary. The top level key is the mode number
        and the bottom level keys are
            * angFreq       - The angular frequency
            * angFreqStd    - The standard deviation of the angular frequency
            * growthRate    - The growth rate
            * growthRateStd - The standard deviation of the growth rate
            * startTime     - Start time of the linear segment
            * endTime       - End time of the linear segment
    """
    #}}}

    # Select modes from 1 to maxMode+1 (exclude the offset mode)
    modes = modes[:, 1:maxMode+1]

    # The absolute value is the magnitude of the signal
    # http://dsp.stackexchange.com/questions/23994/meaning-of-real-and-imaginary-part-of-fourier-transform-of-a-signal
    absModes = np.abs(modes)

    # This is part of the algorithm used to detect the straigth line in
    # the logarithm of the signal
    # Bin the time axis
    binSize = 20
    bins = np.arange(0,len(modes[:,0]),binSize)

    # Finding mode number
    results = {}

    # Transpose in order to loop over the modes
    for modeNr, curMode in enumerate(absModes.transpose()):
        # Place holders for the growth rates and the mean square error
        growthRates  = []
        startIndices = []

        # Loop over the bins, start from 1 as we will index b-1
        for bNr, b in enumerate(range(1, len(bins))):
            # Find the growth rate of the current bin
            startIndex    = bins[b-1]
            endIndex      = bins[b]
            curTime       = time[startIndex: endIndex]
            modeAtCurTime = curMode[startIndex: endIndex]

            growthRate, _ = linRegOfExp(curTime, modeAtCurTime)
            growthRates .append(growthRate)
            startIndices.append(startIndex)

        #{{{Find growth rates from what is defined as straight segments in the plot

        # Initialize the previous growth rate
        prevRate = growthRates[0]

        # The growth rates and the corresponding standar deviation is found by the
        # criteria that the growth rate should stay the same ober at least 4
        # bins. If the criterion is met, a new linear regression will be taken from
        # the max and min of the time indices.

        # Place holders (used in case there are several hits)
        growthRates = []
        sigmaBs     = []
        startTimes  = []
        endTimes    = []
        startIndex = 0
        hits = 0
        for gNr, growthRate in enumerate(growthRates[1:]):
            # If the growth rate is within 1%
            if growthRate*0.99 <= prevRate <= growthRate*1.01:
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
                    endIndex           = bins[int(startIndex/binSize) + hits]
                    curTime            = time[startIndex: endIndex]
                    modeAtCurTime      = curMode[startIndex: endIndex]
                    growthRate, sigmaB = linRegOfExp(curTime, modeAtCurTime)
                    growthRates.append(growthRate)
                    sigmaBs    .append(sigmaB)
                    # We will currently use indices for the start times
                    # as these are easier to deal with when finding the
                    # angular velocity. The indices will be converted to
                    # actual time in the end
                    startTimes.append(startIndex)
                    endTimes  .append(endIndex)
                # Reset the counter
                hits = 0
            prevRate = growthRate

        # If the last element in the growthRates gave a hit
        if hits >= 4:
            endIndex           = bins[int(startIndex/binSize) + hits]
            curTime            = time[startIndex: endIndex]
            modeAtCurTime      = curMode[startIndex: endIndex]
            growthRate, sigmaB = linRegOfExp(curTime, modeAtCurTime)
            growthRates.append(growthRate)
            sigmaBs    .append(sigmaB)
            # We will currently use indices for the start times
            # as these are easier to deal with when finding the
            # angular velocity. The indices will be converted to
            # actual time in the end
            startTimes.append(startIndex)
            endTimes  .append(endIndex)

        if len(growthRates) > 0:
            if len(growthRates) > 1:
                message = ("{0}{1}WARNING: "\
                           "Found {2} different linear segments. "\
                           "Selecting first{1}{0}")
                print(message.format("\n"*2, "!"*5, len(growthRates)))

            # Select the first real growth rate
            results[modeNr]["growthRate"]    = growthRates[0]
            results[modeNr]["growthRateStd"] = sigmaBs    [0]
            results[modeNr]["startIndex"]    = startTimes [0]
            results[modeNr]["endIndex"]      = endTimes   [0]
        else:
            results[modeNr] = None
        #}}}

        #{{{Finding the real part
        # If no growth rate was found
        if results[modeNr] is None:
            continue
        # Remember: The startTime and endTime are currently indices
        startIndex = results[modeNr]["startTime"]
        endIndex   = results[modeNr]["endTime"]
        deltaT     = time[1] - time[0]
        # Create the place holder for the angular frequency
        angularFreq = np.zeros(endIndex-startIndex)
        # Loop over the indices between start and stop
        # startIndex+1 as we will index i-1
        # endIndex+1 to include the last point
        for nr, i in enumerate(range(startIndex+1, endIndex+1)):
            # atan2 in [-pi, pi]
            prevPhaseShift = np.arctan2(curMode[i-1].imag, curMode[i-1].real)
            curPhaseShift  = np.arctan2(curMode[i  ].imag, curMode[i  ].real)
            # phaseShiftDiff in [0, 2*pi]
            phaseShiftDiff = curPhaseShift - prevPhaseShift
            # Ensure that no wrap around has occured
            if max(prevPhaseShift,curPhaseShift)+abs(phaseShiftDiff) > 2*np.pi:
                # Wrap around occured, adding 2*pi to smallest and
                # recalculating
                if prevPhaseShift > curPhaseShift:
                    curPhaseShift += 2*np.pi
                else:
                    prevPhaseShift += 2*np.pi
                # Recalculate the diff
                phaseShiftDiff = curPhaseShift - prevPhaseShift

            # The angular speed (angular frequency) has units rad/s.
            # Remember that if angularFreq*t = 2*pi the perturbation has
            # revolved one time
            angularFreq[nr] = phaseShiftDiff/deltaT

        # Calculate the mean and the spread (the standard deviation), note
        # that numpys std misses a minus 1 in the denominator, but as N is
        # high, this is negligible
        results[modeNr]["angFreq"]    = angularFreq.mean()
        results[modeNr]["angFreqStd"] = angularFreq.std()
        #}}}

        return results
#}}}

#{{{plotGrowthRates
def plotGrowthRates():
    """
    Plots the growth rates

    Parameters
    ----------
FIXME: Implement
    results : dict
        A dictionary of dictionary. The top level key is the mode number
        and the bottom level keys are
            * angFreq       - The angular frequency
            * angFreqStd    - The standard deviation of the angular frequency
            * growthRate    - The growth rate
            * growthRateStd - The standard deviation of the growth rate
            * startTime     - Start time of the linear segment
            * endTime       - End time of the linear segment
    """

    raise NotImplementedError("plotGrowthRates not implemented")
#}}}
