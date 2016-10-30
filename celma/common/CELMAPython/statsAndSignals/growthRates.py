#!/usr/bin/env python

"""
Contains functions to calculate and plot the angular frequency and
growth rate of a time trace of a spatial FFT
"""

from ..plotHelpers import PlotHelper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import os

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
                results[mNr]["allGrowthRates"] = tuple(finalGrowthRates)
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

#{{{PlotGrowthRates
class PlotGrowthRates(object):
    """Class which contains the growth rates and the plotting configuration."""

    #{{{__init___
    def __init__(self                      ,\
                 paths                     ,\
                 growthRates               ,\
                 convertToPhysical = False ,\
                 showPlot          = False ,\
                 savePlot          = False ,\
                 extension         = "png" ,\
                 savePath          = "."   ,\
                 pltSize           = (9,12),\
                 ):
        #{{{docstring
        """
        The constructor for the PlotGrowthRates object.

        Sets the member data.

        Parameters
        ----------
        paths : str
            Paths to collect from (used to make the PlotHelper object)
        growthRates : DataFrame
            The data frame to be plotted. The data frame must be
            orgnized such that:
            1. The x axis will be given by splitting growthRates.columns
               by "=" and take the value right og the "="
            2. The index only has two levels
            3. The plot name is given at the first level of
               growthRates.index
            4. The growth rate, angular frequency and their standar
               deviation is given in the innermost level of growthRates
        showPlot : bool
            If the plots should be displayed.
        savePlot : bool
            If the plots should be saved.
        extension : str
            Extension to use on the plots
        savePath : str
            Path to save destination. Must exist.
        pltSize : tuple
            Size of the plots given as (x, y)
        """
        #}}}

        # Set the member data
        self._df        = growthRates
        self._showPlot  = showPlot
        self._savePlot  = savePlot
        self._extension = extension
        self._savePath  = savePath

        # Make the PlotHelper object
        self._helper = PlotHelper(paths[0][0]                          ,\
                                  useSpatial        = False            ,\
                                  convertToPhysical = convertToPhysical,\
                                 )

        # Set the plot size
        self._pltSize = pltSize

        # Dict used as an equivalent to C++ maps
        self._mapToPltText = {"modeNr":"$\mathrm{Mode number}$",\
                              "B0"    :"$B_0$"                 ,\
                              "Te0"   :"$T_e$"                 ,\
                              "nn"    :"$n_n$"                 ,\
                             }

        self._errorbarOptions = {"color"     :"k",\
                                 "fmt"       :"o",\
                                 "fmt"       :"o",\
                                 "markersize":10 ,\
                                 "ecolor"    :"k",\
                                 "capsize"   :7  ,\
                                 "elinewidth":3  ,\
                                 }

        self._markeredgewidth = 3
        self._all = slice(None)
    #}}}

    #{{{plotGrowthRates
    def plotGrowthRates(self):
        #{{{docstring
        """
        Plots the growth rates.
        """
        #}}}

        # Loop through the figures
        for plotLabel in self._df.index.levels[0]:
            fig = plt.figure(figsize = self._pltSize)
            gs  = GridSpec(nrows=2, ncols=1)

            imAx   = fig.add_subplot(gs[0])
            realAx = fig.add_subplot(gs[1], sharex=imAx)

            plotLabelSplit = plotLabel.split("=")

            # We can now access unsing the loc method
            # The syntax is
            # self._df.loc[("indexLevel1", "indexLevel2", ...),\
            #              ("columnsLevel1", "columnsLevel2", ...)]
            # http://pandas.pydata.org/pandas-docs/stable/cookbook.html#cookbook-selection
            try:
                xAxis = tuple(float(txt.split("=")[1])\
                        for txt in\
                        self._df.loc[(plotLabel, "growthRate"), self._all].\
                        index)
            except KeyError as ke:
                message="{0}{1}WARNING: Only NaNs found in {2}. Skipping{1}{0}"
                print(message.format("\n", "!"*4, ke.args[0]))
                continue

            yAxisIm=self._df.loc[(plotLabel,"growthRate"), self._all].values
            yErrIm =self._df.loc[(plotLabel,"growthRateStd"), self._all].values
            yAxisRe=self._df.loc[(plotLabel,"angFreq"), self._all].values
            yErrRe =self._df.loc[(plotLabel,"angFreqStd"), self._all].values

            indexTxt =\
                self._df.loc[(plotLabel, "growthRate"), self._all].index[0].\
                                                                split("=")[0]

            # Convert units
            plotLabelSplit[-1], plotLabelNorm, plotLabelUnits =\
                self._helper.physicalUnitsConverter(plotLabelSplit[-1],\
                                                    plotLabelSplit[0])
            # NOTE: The whole dataframe can be multiplied using
            #       df.multiply, but we here only multiply the array for
            #       cleaner "helper" code
            xAxis, xAxisNorm, xAxisUnits =\
                self._helper.physicalUnitsConverter(xAxis  , indexTxt)
            yAxisIm, yAxisNorm, yAxisUnits =\
                self._helper.physicalUnitsConverter(yAxisIm, "growthRate")
            yErrIm, _, _  =\
                self._helper.physicalUnitsConverter(yErrIm , "growthRate")
            yAxisRe, _, _ =\
                self._helper.physicalUnitsConverter(yAxisRe, "growthRate")
            yErrRe, _, _  =\
                self._helper.physicalUnitsConverter(yErrRe , "growthRate")

            # Plot the growth rates
            (_, caps, _) = imAx.errorbar(xAxis,\
                                         yAxisIm,\
                                         yerr=yErrIm,\
                                         **self._errorbarOptions)

            # Set capsize
            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            # Angular frequencies
            (_, caps, _) = realAx.errorbar(xAxis,\
                                           yAxisRe,\
                                           yerr=yErrRe,\
                                           **self._errorbarOptions)

            # Set capsize
            for cap in caps:
                cap.set_markeredgewidth(self._markeredgewidth)

            # Add 10% margins for readability
            imAx  .margins(x=0.1, y=0.1)
            realAx.margins(x=0.1, y=0.1)

            rot = None if indexTxt == "modeNr" else 45
            PlotHelper.makePlotPretty(imAx,   yprune = "both", rotation = rot)
            PlotHelper.makePlotPretty(realAx, yprune = "both", rotation = rot)

            # Set the text
            if self._helper.convertToPhysical:
                if plotLabelSplit[0] == "modeNr":
                    suptitle = "{}$={}$".\
                        format(self._mapToPltText[plotLabelSplit[0]],\
                               plotLabelSplit[1])
                else:
                    suptitle = "{}$={}$ $[{}]$".\
                        format(self._mapToPltText[plotLabelSplit[0]],\
                               plotLabelSplit[1],\
                               plotLabelUnits)
                imLabel = r"$\omega_I$ $[{}]$".format(yAxisUnits)
                reLabel = r"$\omega_R$ $[{}]$".format(yAxisUnits)
            else:
                suptitle = "{}$={}{}$".\
                    format(self._mapToPltText[plotLabelSplit[0]],\
                           plotLabelSplit[1],\
                           plotLabelNorm)
                imLabel = r"$\omega_I{}$".format(yAxisNorm)
                reLabel = r"$\omega_R{}$".format(yAxisNorm)

            fig.suptitle(suptitle)
            imAx  .set_ylabel(imLabel)
            realAx.set_ylabel(reLabel)

            imAx.tick_params(labelbottom="off")

            # Set xlabel
            if indexTxt == "modeNr":
                xlabel = self._mapToPltText[indexTxt]
            else:
                if self._helper.convertToPhysical:
                    xlabel = "{} $[{}]$".\
                            format(self._mapToPltText[indexTxt], xAxisUnits)
                else:
                    xlabel = "{}{}".\
                            format(self._mapToPltText[indexTxt], xAxisNorm)

            realAx.set_xlabel(xlabel)

            # Adjust the subplots
            fig.subplots_adjust(hspace=0)

            # Sort the xAxis and yAxis, remove NaN's, and use them as ticks
            xAxis, yAxisIm = zip(*sorted(zip(xAxis, yAxisIm)))
            nonNan = np.where(np.isfinite(yAxisIm))[0]
            # Cut the NaN values from the xAxis (plus 1 as slice
            # excludes the last)
            ticks = tuple(xAxis[nonNan[0]:nonNan[-1]+1])

            # Set the ticks
            realAx.xaxis.set_ticks(ticks)
            imAx  .xaxis.set_ticks(ticks)

            if self._savePlot:
                fileName = "{}.{}".\
                    format(os.path.join(self._savePath,\
                           "growthrate_{}_{}".format(plotLabel, indexTxt)),\
                           self._extension)
                self._helper.savePlot(fig, fileName)

            if self._showPlot:
                plt.show()
    #}}}
#}}}
