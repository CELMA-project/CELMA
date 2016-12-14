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
    """
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
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is used.
        If mode is "fluct" the fluctuations are used.
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

    # Call the time calcTimeTrace function in order to get a time trace
    timeTraces =\
        calcTimeTrace(paths                                ,\
                      varName                              ,\
                      xInd                                 ,\
                      yInd                                 ,\
                      zInd                                 ,\
                      convertToPhysical = convertToPhysical,\
                      mode              = mode             ,\
                      tSlice            = tSlice           ,\
                      )

    YOU ARE HERE!
                    nu      <nu>                        tilde n tilde u
    user choose, normal,  pol avg (collect profile), fluct
    these subtracted gives <<n><u>>
    # To lowest order ExB is the only radial advection
    radialExB = calcTimeTrace(paths                        ,\
                      varName                              ,\
                      xInd                                 ,\
                      yInd                                 ,\
                      zInd                                 ,\
                      convertToPhysical = convertToPhysical,\
                      mode              = mode             ,\
                      tSlice            = tSlice           ,\
                      )


    # Initialize the output
    PSD = {}

    # Make the keys
    xKey = "{}PSDX".format(varName)
    yKey = "{}PSDY".format(varName)

    # Obtain the PSD
    for key in timeTraces.keys():
        # Initialize the PSD
        PSD[key] = {}

        # Sampling frequency
        fs = 1/(timeTraces[key]["time"][1] - timeTraces[key]["time"][0])

        # window = None => window = "boxcar"
        # scaling = density gives the correct units
        PSD[key][xKey], PSD[key][yKey] =\
            periodogram(timeTraces[key],\
                        fs=fs, window=None, scaling="density")

    # NOTE: If timeTraces was converted to physical units, then PSD is
    #       in physical units as well
    return PSD
#}}}
