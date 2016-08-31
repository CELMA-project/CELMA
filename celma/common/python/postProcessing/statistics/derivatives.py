#!/usr/bin/env python

"""
Contains derivative functions
"""

import numpy as np
from scipy.fftpack import diff
from boutdata import collect

#{{{DDZ1D
def DDZ1D(poloidalProfile, xInd, yInd):
    """
    Calculates the first theta-derivative of a profile using spectral
    methods (assuming that we are using cylinder geometry).

    Note that in cylinder geometry, this derivative has to be multiplied
    with 1/J. Also note that the profile must go from [0, 2*pi[ (which
    is standard output from BOUT++, and NOT [0, 2*pi]

    Parameters
    ----------
    poloidalProfile : array
        A 1-D array of a profile (fixed xInd [rho index] and yInd
        [cylinder length index])

    xInd: int
        What x-index to use in J
    yInd: int
        What y-index to use in the dx

    Returns
    -------
    out : iterable
        The theta derivative of the profile
    """


    J = collect("J", xind = [xInd, xInd],
                     yind = [yInd, yInd],
                     xguards=False, yguards=False)

    out = diff(poloidalProfile)*J

    return out
#}}}

#{{{DDX1D
def DDX1D(profile, yInd):
    """
    Calculates the first rho-derivative of a profile using a second order
    stencil (assuming that we are using cylinder geometry).

    Parameters
    ----------
    profile : array
        A 1-D array of a profile (fixed yInd [cylinder length index] and zInd
        [theta index])
    yInd: int
        What y-index to use in the dx

    Returns
    -------
    out : iterable
        The rho derivative of the profile
    """

    dx = collect("dx", xguards=True, yguards=False)[:, yInd]

    out = np.zeros(profile.shape)

    out[0] = out[-1] = np.NaN
    # Pad with two ghost points
    for xInd in range(1, len(profile)-1):
        # 2nd order scheme
        out[xInd] = 0.5*(-profile[xInd-1] + profile[xInd+1])/dx[xInd, yInd]

    return out
#}}}

#{{{getRadialExBVelocityZPlane
def getRadialExBVelocityZPlane(phi, yInd):
    """
    Returns the radial velocity at a certain z-plane (i.e. a plane perpendicular
    to B)

    The ExB velocity for a Clebsch system can be found in section B.5.1
    in the BOUT++ coordinates manual. However, the cylindrical
    coordinate system is not a Clebsch system, but the metric overlaps.
    In order to get the cylindrical coordinates ExB, we must multiply
    the ExB velocity in BOUT++ with B (i.e. divide by rho). Thus, the
    radial ExB velocity is the cylindrical theta derivative of phi.

    Parameters
    ----------
    phi : array
        The potential at a certain z-plane (i.e. for a specified yInd)
    yInd : int
        The yInd for the phi array

    Returns
    -------
    radialExB : array
        The ExB velocity a certain z-plane (i.e. for a specified yInd)
    """

    if len(phi.shape) == 3:
        # The time dimension is present
        nXinds = phi.shape[1]
    elif len(phi.shape) == 2:
        # The time dimension is not present
        nXinds = phi.shape[0]

    radialExB = np.zeros(phi.shape)

    for xInd in range(nXinds):
        radialExB[xInd, :] = DDZ1D(phi, xInd, yInd)

    return radialExB
#}}}

#{{{findLargestGrad
def findLargestGrad(var, yInd = 0):
    """
    Returns the value of the maximum gradient and the corresponding
    index

    Parameters
    ----------
    var : array
        The variable to investigate.
    yInd : int
        The index to evaluate dx in

    Returns
    -------
    maxGrad : float
        The value of the maximum gradient
    maxGradInd : int
        The index at which the maximum gradient can be found at
    """

    MXG = collect("MXG", xguards=False, yguards=False)

    ddxVar = DDX1D(var, yInd)
    # Subtract MXG as xguards are collected
    maxGradInd = int(np.where(np.abs(ddxVar) == np.max(np.abs(ddxVar)))[0]) - MXG

    maxGrad = ddxVar[maxGradInd]

    return maxGrad, maxGradInd
#}}}
