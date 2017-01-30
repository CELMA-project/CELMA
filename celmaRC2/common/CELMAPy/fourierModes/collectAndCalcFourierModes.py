#!/usr/bin/env python

"""
Contains class for collecting and calculating the fourier modes
"""

from ..superClasses import CollectAndCalcPointsSuperClass
from ..collectAndCalcHelpers import (collectTime,\
                                     collectPoloidalProfile,\
                                     calcN,\
                                     calcUIPar,\
                                     calcUEPar,\
                                     slicesToIndices)
import numpy as np

#{{{CollectAndCalcFourierModes
class CollectAndCalcFourierModes(CollectAndCalcPointsSuperClass):
    """
    Class for collecting and calcuating the fourier modes
    """

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates the fourier modess.

        Returns
        -------
        fourierModes : dict
            Dictionary where the keys are on the form "rho,z".
            The value is a dict containing of
            {varName:fourierModes, "time":time}.
            The fourierModes is a 4d array.
        """
        #}}}

        # Guard
        if len(self._notCalled) > 0:
            message = "The following functions were not called:\n{}".\
                        format("\n".join(self._notCalled))
            raise RuntimeError(message)

        # Initialize output
        fourierModes = {}
        tCounter = 0

        for x, y in zip(self._xInd, self._yInd):
            # NOTE: The indices
            rho = self._dh.rho[x]
            par = self._dh.z  [y]

            # Add key and dict to fourierModes
            key = "{},{}".format(rho,par)
            fourierModes[key] = {}

            if self._tSlice is not None:
                t = slicesToIndices(self._collectPaths, self._tSlice[tCounter], "t")
                tStep = self._tSlice[tCounter].step
            else:
                t = None
                tStep = None
            tCounter += 1

            var, time = self._collectWrapper(fourierModes,key,x,y,t,tStep)

            if self.uc.convertToPhysical:
                fourierModes[key][self._varName] =\
                        self.uc.physicalConversion(var , self._varName)
                fourierModes[key]["time"]  =\
                        self.uc.physicalConversion(time, "t")
            else:
                fourierModes[key][self._varName] = var
                fourierModes[key]["time"]        = time

        return fourierModes
    #}}}

    #{{{convertTo2D
    def convertTo2D(self, fourierModes):
        #{{{docstring
        """
        Converts the 4d array to a 2d array.

        Parameters
        ----------
        fourierModes : dict
            Output from executeCollectAndCalc.

        Returns
        -------
        fourierModes : dict
            As the input, but 2d traces on the form (t, nz) rather than 4d.
        """
        #}}}

        for key in fourierModes.keys():
            fourierModes[key][self._varName] =\
                    fourierModes[key][self._varName][:,0,0,:]

        return fourierModes
    #}}}

    #{{{_collectWrapper
    def _collectWrapper(self,fourierModes,key,x,y,t,tStep):
        #{{{docstring
        """
        Collects the variable and the time.

        If the varName is n, uIPar or uEPar, calculation will be done
        through _calcNonSolvedVars

        Parameters
        ----------
        fourierModes : dict
            Dict containing the fourier modess
        key : str
            Key with the point position
        x : int
            The x index to collect from
        y : int
            The y index to collect from
        z : int
            The z index to collect from
        t : [None|tuple]
            The collect-like slice in t
        tStep : [None|int]
            The step to chop the time variable in

        Returns
        -------
        var : array-4d
            The fourier transformed collected array.
        time : array
            The time array.
        """
        #}}}

        time = collectTime(self._collectPaths, tInd=t)
        collector = collectPoloidalProfile
        indArgs = (x, y)
        if not(self._varName == "n" or\
               self._varName == "uIPar" or\
               self._varName == "uEPar"):
            var = collector(self._collectPaths,\
                            self._varName,\
                            *indArgs, tInd=t)
        else:
            var = self._calcNonSolvedVars(collector,indArgs,t)

        # Fourier transform
        var = np.fft.fft(var)

        if tStep is not None:
            # Slice the variables with the step
            # Make a new slice as the collect dealt with the start and
            # the stop of the slice
            newSlice = slice(None, None, tStep)

            var  = var [newSlice]
            time = time[newSlice]

        return var, time
    #}}}

    #{{{_calcNonSolvedVars
    def _calcNonSolvedVars(self,collector,indArgs,t):
        #{{{docstring
        """
        Calculates variables wich are not solved in the simulation

        NOTE: We set normalized True here as the collected
              variables are normalized.
              Conversion to physical happens later in
              executeCollectAndCalc.

        Parameters
        ----------
        collector : func
            Function to use for collection
        indArgs : tuple
            Tuple of index arguments.
            Either x and y.

        Returns
        -------
        var : 4d-array
            The fourier modes of the calculated variable.
        """
        #}}}

        normalized = True

        lnN = collector(self._collectPaths,\
                        "lnN"             ,\
                        *indArgs, tInd=t)
        n = calcN(lnN, normalized, uc = self.uc)
        if self._varName == "n":
            return n
        else:
            momDensPar = collector(self._collectPaths,\
                                   "momDensPar"      ,\
                                   *indArgs, tInd=t)
            uIPar = calcUIPar(momDensPar, n)
            if self._varName == "uIPar":
                return uIPar
            else:
                jPar = collector(self._collectPaths,\
                                "jPar"            ,\
                                *indArgs, tInd=t)
                uEPar = calcUEPar(uIPar     ,\
                                  jPar      ,\
                                  n         ,\
                                  normalized)
                return uEPar
    #}}}

    #{{{calcMagnitude
    def calcMagnitude(self, fourierModes2d):
        #{{{docstring
        """
        Calculates the magnitude of the 2d fourier signal.

        Parameters
        ----------
        fourierModes2d : dict
            Output from executeCollectAndCalc and _convertTo2D.

        Returns
        -------
        fourierModes2d : dict
            As the input, but contains the key varNameMagnitude with
            values on the form (t, (nz/2) + 1) for each position.
        """
        #}}}

        for key in fourierModes2d.keys():
            modes = fourierModes2d[key][self._varName].copy()
            tSize, N  = modes.shape
            nyquistMode = int(N/2) + 1
            magnitude = np.zeros((tSize, nyquistMode))
            #{{{ NOTE: We are dealing with a real signal:
            #          As the fourier transform breaks the signal up in
            #          cisoids there will be one part of the signal in
            #          the positive rotating ciscoid and one in the
            #          negative (negative frequencies) for a given mode
            #          number. We need to take into account both in
            #          order to calculate the amplitude.
            # http://dsp.stackexchange.com/questions/431/what-is-the-physical-significance-of-negative-frequencies?noredirect=1&lq=1
            # http://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
            #}}}
            # Magnitude of the signal
            # https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Definition
            # The offset mode and the Nyquist mode
            magnitude[:,0]             = np.abs(modes[:,0])/N
            # Minus 1 as indices count from 0
            magnitude[:,nyquistMode-1] = np.abs(modes[:,nyquistMode])/N
            # Loop over all the modes
            # NOTE: Range exludes the last point
            # Minus 1 as indices count from 0
            for modeNr in range(1, nyquistMode-1):
                posFreq              = np.abs(modes[:,  modeNr])
                negFreq              = np.abs(modes[:, -modeNr])
                magnitude[:, modeNr] = (posFreq + negFreq)/N

            # Insert into the dict
            fourierModes2d[key][self._varName+"Magnitude"] = magnitude
        return fourierModes2d
    #}}}

    #{{{calcAngularFrequency
    def calcAngularFrequency(self, fourierModes2d):
        #{{{docstring
        """
        Calculates the phaseShift of the 2d fourier signal.

        Parameters
        ----------
        fourierModes2d : dict
            Output from executeCollectAndCalc and _convertTo2D.

        Returns
        -------
        fourierModes2d : dict
            As the input, but contains the key varNameAngularFrequency
            with values on the form (t, (nz/2) + 1) for each position.
            NOTE: A negative angular frequency means that the
                  perturbations are moving in the negative theta
                  direction.
        """
        #}}}

        for key in fourierModes2d.keys():
            modes    = fourierModes2d[key][self._varName].copy()
            time     = fourierModes2d[key]["time"]
            tSize, N = modes.shape
            nyquistMode = int(N/2) + 1
            angularFreq = np.zeros((tSize-1, nyquistMode))
            #{{{ NOTE: We are dealing with a real signal:
            #          As the signal is real only one of the phase sifts
            #          are needed. Notice that for a real signal the
            #          imaginary part occurs as a complex conjugate pair
            # http://dsp.stackexchange.com/questions/431/what-is-the-physical-significance-of-negative-frequencies?noredirect=1&lq=1
            # http://dsp.stackexchange.com/questions/4825/why-is-the-fft-mirrored
            #}}}
            # The phase shift is found from atan2
            # http://dsp.stackexchange.com/questions/23994/meaning-of-real-and-imaginary-part-of-fourier-transform-of-a-signal
            # atan2 in [-pi, pi]
            for modeNr in range(0, nyquistMode):
                for tInd in range(1, tSize):
                    deltaT = time[tInd] - time[tInd-1]
                    prevPhaseShift = np.arctan2(modes[tInd-1, modeNr].imag,\
                                                modes[tInd-1, modeNr].real)
                    curPhaseShift  = np.arctan2(modes[tInd  , modeNr].imag,\
                                                modes[tInd  , modeNr].real)
                    # phaseShiftDiff in [0, 2*pi]
                    phaseShiftDiff = prevPhaseShift - curPhaseShift

                    # Corrections
                    if curPhaseShift*prevPhaseShift < 0\
                       and abs(curPhaseShift) + abs(prevPhaseShift) > np.pi:
                        if curPhaseShift < 0:
                            # We are going from pi to -pi
                            # In order to avoid the discontinuity, we turn
                            # curPhaseShift and prevPhaseShift to the opposite
                            # quadrants
                            tempCurPhase   = np.pi + curPhaseShift
                            tempPrevPhase  = np.pi - prevPhaseShift
                            phaseShiftDiff = -(tempCurPhase + tempPrevPhase)
                        else:
                            # We are going from -pi to pi
                            # In order to avoid the discontinuity, we
                            # turn curPhaseShift and prevPhaseShift to
                            # the opposite quadrants
                            tempCurPhase   = np.pi - curPhaseShift
                            tempPrevPhase  = np.pi + prevPhaseShift
                            phaseShiftDiff = tempCurPhase + tempPrevPhase

                    # The angular speed (angular frequency) has units rad/s.
                    # Remember that if angularFreq*t = 2*pi the perturbation has
                    # revolved one time
                    angularFreq[tInd-1, modeNr] = phaseShiftDiff/deltaT

            # Insert into the dict
            fourierModes2d[key][self._varName+"AngularFrequency"] = angularFreq
        return fourierModes2d
    #}}}
#}}}
