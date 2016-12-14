#!/usr/bin/env python

"""
Contains the calcTimeTrace calculation
"""

from .averages import polAvg
from ..plotHelpers import (PlotHelper,\
                           collectiveCollect,\
                           safeCollect,\
                           DDZ,\
                           findLargestRadialGrad)
import numpy as np
from scipy.stats import kurtosis, skew
from scipy.signal import periodogram

#{{{calcTimeTrace
def calcTimeTrace(paths                      ,\
                  varName                    ,\
                  xInd                       ,\
                  yInd                       ,\
                  zInd                       ,\
                  convertToPhysical = True   ,\
                  mode              = "fluct",\
                  tSlice            = None):
    #{{{docstring
    """
    Function which calculates the time traces

    Parameters
    ----------
    paths : tuple of strings
        The paths to collect from
    varName : str
        Variable to collect
    xInd : tuple of ints
        A tuple of the xInds to collect use when collecting.
    yInd : tuple of ints
        The same as xInd, but for the y-index.
    zInd : tuple of ints
        The same as xInd, but for the z-index.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    tSlice : [None|Slice}
        Whether or not to slice the time trace

    Returns
    -------
    timeTraces : dict
        Dictionary where the keys are on the form "rho,theta,z".
        The value is a dict containing of
        {varName:timeTrace, "time":time}
    """
    #}}}

    # Create the units convertor object
    uc = UnitsConverter(paths[0], convertToPhysical)
    # Toggle convertToPhysical in case of errors
    convertToPhysical = uc.convertToPhysical
    # Create the dimensions helper object
    self.dh = DimensionsHelper(paths[0], uc)

    timeTraces = {}
    tCounter = 0
    for x, y, z in zip(xInd, yInd, zInd):
        # NOTE: The indices
        rho   = dh.rho  (int(x))
        theta = dh.theta(int(z))
        z     = dh.z    (int(y))

        # Add key and dict to timeTraces
        key = "{},{},{}".format(x,y,z)
        timeTraces[key] = {}

        if tSlice is not None:
            tStart = tSlice[tCounter].start
            tEnd   = tSlice[tCounter].end
            tCounter += 1
        else:
            tStart = None
            tEnd   = None

        t = (tStart, tEnd)

        var, time = collectPointTime(paths, varName, x, y, z, t)

        # Reshape
        var = var.flatten()

        if tSlice is not None:
            # Slice the variables with the step
            # Make a new slice as the collect dealt with the start and
            # the stop of the slice
            newSlice = slice(None, None, tSlice.step)

            var  = var [newSlice]
            time = time[newSlice]

        if mode == "fluct"
            var = var - var.mean()
        elif mode != "normal":
            raise NotImplementedError("'{}'-mode not implemented".format(mode))

        if convertToPhysical:
            timeTraces[key][varName] = uc.physicalConversion(var , varName)
            timeTraces[key]["time"]  = uc.physicalConversion(time, "t")
        else:
            timeTraces[key][varName] = var
            timeTraces[key]["time"]  = time

    return timeTraces
#}}}
