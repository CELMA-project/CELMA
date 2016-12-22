#!/usr/bin/env python

"""
Contains the calcTimeTrace calculation
"""

from ..calcHelpers import (polAvg,\
                           collectPoint,\
                           collectTime,\
                           collectPoloidalProfile,\
                           DimensionsHelper)
from ..unitsConverter import UnitsConverter


#{{{calcTimeTrace
def calcTimeTrace(*args, **kwargs):
    #{{{docstring
    """
    Wrapper function for calcTimeTrace4d.

    This function will which recasts to 1d, and pop eventual zInd.

    Parameters
    ----------
    See calcTimeTrace4d for details.

    Returns
    -------
    See calcTimeTrace4d for details.
    """
    #}}}

    timeTraces = calcTimeTrace4d(*args, **kwargs)

    varName = args[1]
    mode    = kwargs["mode"]

    for key in timeTraces.keys():
        if mode == "fluct":
            # The fluctuations does not have a specified z
            z = timeTraces[key].pop("zInd")
            timeTraces[key][varName] = timeTraces[key][varName][:,:,:,z:z+1]

        # Reshape
        timeTraces[key][varName] = timeTraces[key][varName].flatten()

    return timeTraces
#}}}

#{{{calcTimeTrace4d
def calcTimeTrace4d(paths                      ,\
                    varName                    ,\
                    xInd                       ,\
                    yInd                       ,\
                    zInd                       ,\
                    convertToPhysical = True   ,\
                    mode              = "fluct",\
                    tSlice            = None   ,\
                    uc                = None   ,\
                    dh                = None   ,\
                    ):
    #{{{docstring
    """
    Function which calculates the time traces, and gives a 4d output

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
    tSlice : [None|Slice]
        Whether or not to slice the time trace
    uc : [None|UnitsConverter]
        If not given, the function will create the instance itself.
        However, there is a possibility to supply this in order to
        reduce overhead.
    dh : [None|DimensionsHelper]
        If not given, the function will create the instance itself.
        However, there is a possibility to supply this in order to
        reduce overhead.

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

    if uc is None:
        # Create the units convertor object
        uc = UnitsConverter(paths[0], convertToPhysical)
    # Toggle convertToPhysical in case of errors
    convertToPhysical = uc.convertToPhysical

    if dh is None:
        # Create the dimensions helper object
        dh = DimensionsHelper(paths[0], uc)

    timeTraces = {}
    tCounter = 0
    for x, y, z in zip(xInd, yInd, zInd):
        # NOTE: The indices
        rho   = dh.rho      [x]
        theta = dh.thetaDeg[z]
        par   = dh.z        [y]

        # Add key and dict to timeTraces
        key = "{},{},{}".format(rho,theta,par)
        timeTraces[key] = {}

        if tSlice is not None:
            tStart = tSlice[tCounter].start
            tEnd   = tSlice[tCounter].end
            t = (tStart, tEnd)
        else:
            t = None

        tCounter += 1

        time = collectTime(paths, tInd=t)

        if mode == "normal":
            var = collectPoint(paths, varName, x, y, z, tInd=t)
        elif mode == "fluct":
            var = collectPoloidalProfile(paths, varName, x, y, tInd=t)
            var = (var - polAvg(var))
            timeTraces[key]["zInd"] = z
        else:
            raise NotImplementedError("'{}'-mode not implemented".format(mode))

        if tSlice is not None:
            # Slice the variables with the step
            # Make a new slice as the collect dealt with the start and
            # the stop of the slice
            newSlice = slice(None, None, tSlice.step)

            var  = var [newSlice]
            time = time[newSlice]

        if convertToPhysical:
            timeTraces[key][varName] = uc.physicalConversion(var , varName)
            timeTraces[key]["time"]  = uc.physicalConversion(time, "t")
        else:
            timeTraces[key][varName] = var
            timeTraces[key]["time"]  = time

    return timeTraces
#}}}
