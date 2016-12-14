#!/usr/bin/env python

"""
Contains the radial flux calculation
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

#{{{calcRadialFlux
def calcRadialFlux(paths                      ,\
                   varName                    ,\
                   xInd                       ,\
                   yInd                       ,\
                   zInd                       ,\
                   convertToPhysical = True   ,\
                   mode              = "fluct",\
                   tSlice            = None):
    #{{{docstring
    r"""
    Function which calculates the radial flux.

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
    mode : ["normal"|"avg"|"fluct"]
        If mode is "normal" the output is on the form nu.
        If mode is "avg" the output is on the form <nu>.
        If mode is "fluct" the output is on the form \tilde{n}\tilde{u}.
        Note that
        <nu> = <(<n>+\tidle{n})(<u>+\tidle{u})>
             = <<n><u>> + <\tidle{n}\tidle{u})>
             = <n><u> + <\tidle{n}\tidle{u})>
        So that <n><u> is given by the "avg" - "fluct"
    tSlice : [None|Slice}
        Whether or not to slice the time trace

    Returns
    -------
    radialFLux : dict
        Dictionary where the keys are on the form "rho,theta,z".
        The value is a dict containing of
        {varRadialFlux:radialFlux, "time":time}
    """
    #}}}

    uc = UnitsConverter(paths[0], convertToPhysical)
    dh = DimensionsHelper(paths[0], uc)

    # Set mode to call the functions with
    if mode == "avg":
        callMode = "normal"
    else:
        callMode = mode

    # Call the time calcTimeTrace function in order to get a time trace
    varTimeTraces =\
        calcTimeTrace4d(paths                                ,\
                        varName                              ,\
                        xInd                                 ,\
                        yInd                                 ,\
                        zInd                                 ,\
                        convertToPhysical = convertToPhysical,\
                        mode              = callMode         ,\
                        tSlice            = tSlice           ,\
                        uc                = uc               ,\
                        dh                = dh               ,\
                        )

    # To lowest order ExB is the only radial advection
    radialExB =\
        calcTimeTraceRadialDerivative(\
            paths                                ,\
            "phi"                                ,\
            xInd                                 ,\
            yInd                                 ,\
            zInd                                 ,\
            convertToPhysical = convertToPhysical,\
            mode              = callMode         ,\
            tSlice            = tSlice           ,\
            uc                = uc               ,\
            dh                = dh               ,\
            )

    # Initialize the output
    radialFlux = {}
    varRadialFlux = "{}RadialFlux".format(varName)

    keys = sort(radialExB.keys())

    for key in keys():
        raidalFlux[key] = {}
        if mode == "normal":
            radialFlux[key][varRadialFlux] =\
                radialExB[key]*varTimeTraces[key][varName]
        elif mode == "avg":
            radialFlux[key][varRadialFlux] =\
                polAvg(radialExB[key]*varTimeTraces[key][varName])
        elif mode == "fluct":
            fluctExB = radialExB[key] - polAvg(radialExB[key])
            fluctVar = varTimeTraces[key][varName] -\
                       polAvg(varTimeTraces[key][varName])

            z = varTimeTraces[key]["zInd"].pop()
            radialFlux[key][varRadialFlux] = (fluctExB*fluctVar)[:,:,:,z:z]
        else:
            raise NotImplementedError("'{}'-mode not implemented")

        # Flatten
        radialFlux[key][varRadialFlux].flatten()
        radialFlux[key]["time"] = varTimeTraces[key]["time"]

    # NOTE: If varTimeTraces and radialExB was converted to physical units,
    #       then radialFlux is as well
    return radialFlux
#}}}

#{{{calcTimeTraceRadialDerivative
def calcTimeTraceRadialDerivative(paths                      ,\
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
    DDZVar : dict
        Dictionary where the keys are on the form "rho,theta,z".
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

    DDZVars  = {}
    tCounter = 0
    for x, y, z in zip(xInd, yInd, zInd):
        # NOTE: The indices
        rho   = dh.rho  [x]
        theta = dh.theta[z]
        z     = dh.z    [y]

        # Add key and dict to timeTraces
        key = "{},{},{}".format(x,y,z)
        timeTraces[key] = {}

        if tSlice is not None:
            tStart = tSlice[tCounter].start
            tEnd   = tSlice[tCounter].end
        else:
            tStart = None
            tEnd   = None

        tCounter += 1
        t = (tStart, tEnd)

        if mode == "normal":
            var, _ = collectPoloidalProfileTime(paths, varName, x, y, tInd=t)
            DDZVar = DDZ(var, dh.rho[x])
        elif mode == "fluct":
            var, _ = collectPoloidalProfileTime(paths, varName, x, y, tInd=t)
            DDZVar = DDZ(var, dh.rho[x])
            DDZVar = (DDZVar - polAvg(DDZVar))
        else:
            raise NotImplementedError("'{}'-mode not implemented".format(mode))

        if tSlice is not None:
            # Slice the variables with the step
            # Make a new slice as the collect dealt with the start and
            # the stop of the slice
            newSlice = slice(None, None, tSlice.step)

            DDZVar = DDZVar[newSlice]

        if convertToPhysical:
            # NOTE: 1/J is multiplied with DDZ (which has units of angle)
            #       We therefore convert to physical units by first
            #       mutliplying with the factor of var, then with the
            #       factor of 1/rho
            DDZVar = uc.physicalConversion(DDZVar, varName)
            # To convert to 1/J, we use the normFactor of "length"
            DDZVars[key] = uc.normalizedConversion(DDZVar, "length")

    return DDZVars
#}}}
