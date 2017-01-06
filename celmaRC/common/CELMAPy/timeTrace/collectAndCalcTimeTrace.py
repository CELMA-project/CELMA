#!/usr/bin/env python

"""
Contains class for collecting and calculating the time traces
"""

from ..superClasses import CollectAndCalcPointsSuperClass
from ..collectAndCalcHelpers import (polAvg,\
                                     collectPoint,\
                                     collectTime,\
                                     collectPoloidalProfile,\
                                     calcN,\
                                     calcUIPar,\
                                     calcUEPar)

#{{{CollectAndCalcTimeTrace
class CollectAndCalcTimeTrace(CollectAndCalcPointsSuperClass):
    """
    Class for collecting and calcuating the time traces
    """

    #{{{constructor
    def __init__(self            ,\
                 *args           ,\
                 mode  = "normal",\
                 **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor
            * Set member data

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        mode : ["normal"|"fluct"]
            If mode is "normal" the raw data is given as an output.
            If mode is "fluct" the fluctuations are given as an output.
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}

        # Guard
        implemented = ("normal", "fluct")
        if not(mode in implemented):
            message = "mode '{}' not implemented".format(mode)
            raise NotImplementedError(message)

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._mode = mode
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates the time traces.

        Returns
        -------
        timeTraces : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:timeTrace, "time":time}.
            And additional key "zInd" will be given in addition to varName
            and "time" if mode is set to "fluct".
            The timeTrace is a 4d array.
        """
        #}}}

        # Guard
        if len(self._notCalled) > 0:
            message = "The following functions were not called:\n{}".\
                        format("\n".join(self._notCalled))
            raise RuntimeError(message)

        # Initialize output
        timeTraces = {}
        tCounter = 0
        for x, y, z in zip(self._xInd, self._yInd, self._zInd):
            # NOTE: The indices
            rho   = self._dh.rho     [x]
            theta = self._dh.thetaDeg[z]
            par   = self._dh.z       [y]

            # Add key and dict to timeTraces
            key = "{},{},{}".format(rho,theta,par)
            timeTraces[key] = {}

            t = slicesToIndices(self._collectPaths, self._tSlice[tCounter], "t")
            tCounter += 1

            var, time = self._collectWrapper(timeTraces,key,x,y,z,t)

            if self.uc.convertToPhysical:
                timeTraces[key][self._varName] =\
                        self.uc.physicalConversion(var , self._varName)
                timeTraces[key]["time"]  =\
                        self.uc.physicalConversion(time, "t")
            else:
                timeTraces[key][self._varName] = var
                timeTraces[key]["time"]        = time

        return timeTraces
    #}}}

    #{{{convertTo1D
    def convertTo1D(self, timeTraces):
        #{{{docstring
        """
        Converts the 4d array to a 1d array and pops any eventual zInd.

        Parameters
        ----------
        timeTraces : dict
            Output from executeCollectAndCalc.

        Returns
        -------
        timeTraces : dict
            As the input, but 1d traces rather than 4d, and no zInd.
        """
        #}}}

        for key in timeTraces.keys():
            if self._mode == "fluct":
                # The fluctuations does not have a specified z
                z = timeTraces[key].pop("zInd")
                timeTraces[key][self._varName] =\
                        timeTraces[key][self._varName][:,:,:,z:z+1]

            # Reshape
            timeTraces[key][self._varName] =\
                timeTraces[key][self._varName].flatten()

        return timeTraces
    #}}}

    #{{{_collectWrapper
    def _collectWrapper(self,timeTraces,key,x,y,z,t):
        #{{{docstring
        """
        Collects the variable and the time.

        If the varName is n, uIPar or uEPar, calculation will be done
        through _calcNonSolvedVars

        Parameters
        ----------
        timeTraces : dict
            Dict containing the time traces
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

        Returns
        -------
        var : array-4d
            The collected array.
        time : array
            The time array.
        """
        #}}}

        time = collectTime(self._collectPaths, tInd=t)

        if self._mode == "normal":
            collector = collectPoint
            indArgs = (x, y, z)
        elif self._mode == "fluct":
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

        if self._mode == "fluct":
            var = (var - polAvg(var))
            timeTraces[key]["zInd"] = z

        if self._tSlice is not None:
            # Slice the variables with the step
            # Make a new slice as the collect dealt with the start and
            # the stop of the slice
            newSlice = slice(None, None, self._tSlice.step)

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
            Either x,y and z, or x and y.

        Returns
        -------
        var : 4d-array
            The time trace of the calculated variable.
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

    #{{{getDh
    def getDh(self):
        """
        Returns the dimension helper.

        Returns
        -------
        dh : DimensionHelper
            The dimension helper object
        """
        return self._dh
    #}}}

    #{{{getSlices
    def getSlices(self):
        """
        Returns the slices.

        Returns
        -------
        slices : tuple of tuples
            The slices are on the format
            ((xInd1,...), (yInd1,...), (zInd1,...), (tSlice1,...))
        """
        slices = (self._xInd, self._yInd, self._zInd, self._tSlice)
        return slices
    #}}}
#}}}
