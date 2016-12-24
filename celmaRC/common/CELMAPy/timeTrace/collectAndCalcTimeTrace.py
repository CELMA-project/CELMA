#!/usr/bin/env python

"""
Contains class for collecting and calculating the time traces
"""

from ..superClasses import CollectAndCalcPointsSuperClass
from ..collectAndCalcHelpers import (polAvg,\
                                     collectPoint,\
                                     collectTime,\
                                     collectPoloidalProfile,\
                                     DimensionsHelper)

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
            rho   = dh.rho     [x]
            theta = dh.thetaDeg[z]
            par   = dh.z       [y]

            # Add key and dict to timeTraces
            key = "{},{},{}".format(rho,theta,par)
            timeTraces[key] = {}

            if self._tSlice is not None:
                tStart = tSlice[tCounter].start
                tEnd   = tSlice[tCounter].end
                t = (tStart, tEnd)
            else:
                t = None

            tCounter += 1

            var, time self._collectWrapper(x,y,z,t)

            if self.uc.convertToPhysical:
                timeTraces[key][varName] =\
                        self.uc.physicalConversion(var , varName)
                timeTraces[key]["time"]  =\
                        self.uc.physicalConversion(time, "t")
            else:
                timeTraces[key][varName] = var
                timeTraces[key]["time"]  = time

        return timeTraces
    #}}}

    #{{{convertTo1D
    def convertTo1D(timeTraces):
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
    def _collectWrapper(self,x,y,z,t):
        #{{{docstring
        """
        Collects the variable and the time.

        Parameters
        ----------
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

        time = collectTime(collectPaths, tInd=t)

        if self._mode == "normal":
            var = collectPoint(self._collectPaths,\
                               self._varName,\
                               x, y, z, tInd=t)
        elif self._mode == "fluct":
            var = collectPoloidalProfile(self._collectPaths,\
                                         self._varName,\
                                         x, y, tInd=t)
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
#}}}
