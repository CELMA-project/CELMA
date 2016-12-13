#!/usr/bin/env python

"""
Contains derivative functions
"""

import numpy as np
from scipy.fftpack import diff

#{{{DDX
def DDX(var, dx):
    """
    Calculates the first rho-derivative of a profile using a second order
    stencil (assuming that we are using cylinder geometry).

    Parameters
    ----------
    var : array
        Variable to take the first derivative of (must include ghost
        points)
    dx : array
        The grid spacing in x for the different evaluation points

    Returns
    -------
    out : iterable
        The rho derivative of the profile
    """
    if len(var.shape) != 4:
       raise ValueError("Input variable must be 4-dimensional")

    tLen, _, yLen, zLen = var.shape

    out = np.zeros(var.shape)

    for tInd in range(tLen):
        for yInd in range(yLen):
            for zInd in range(zLen):
                # 2nd order scheme
                out[tInd, :, yInd, zInd] =\
                    np.gradient(out[tInd, :, yInd, zInd], dx, edge_order=2)

    return out
#}}}

#{{{DDY
def DDY(var, dy):
    """
    Calculates the first parallel-derivative of a profile using a second order
    stencil (assuming that we are using cylinder geometry).

    Parameters
    ----------
    var : array
        Variable to take the first derivative of (must include ghost
        points)
    dy : array
        The grid spacing in x for the different evaluation points

    Returns
    -------
    out : iterable
        The rho derivative of the profile
    """
    if len(var.shape) != 4:
       raise ValueError("Input variable must be 4-dimensional")

    tLen, xLen, _, zLen = var.shape

    out = np.zeros(var.shape)

    for tInd in range(tLen):
        for xInd in range(xLen):
            for zInd in range(zLen):
                # 2nd order scheme
                out[tInd, xInd, :, zInd] =\
                    np.gradient(out[tInd, xInd, :, zInd], dy, edge_order=2)

    return out
#}}}

#{{{DDZ
def DDZ(var, J):
    """
    Calculates the first theta-derivative of a profile using spectral
    methods (assuming that we are using cylinder geometry).

    Note that in cylinder geometry, this derivative has to be multiplied
    with 1/J. Also note that the profile must go from [0, 2*pi[ (which
    is standard output from BOUT++, and NOT [0, 2*pi]

    Parameters
    ----------
    var : array
        The variable to take the z-derivative of
    J : array
        The Jacobian

    Returns
    -------
    out : iterable
        The theta derivative of the profile
    """
    if len(var.shape) != 4:
       raise ValueError("Input variable must be 4-dimensional")

    tLen, xLen, yLen, _ = var.shape

    out = np.zeros(var.shape)

    for tInd in range(tLen):
        for xInd in range(xLen):
            for yInd in range(yLen):
                out[tInd, xInd, yInd, :] =\
                    diff(var[tInd, xInd, yInd, :])*J[xInd, yInd]

    return out
#}}}

#{{{findLargestRadialGrad
def findLargestRadialGrad(var, dx, MXG = None):
    """
    Returns the value of the maximum gradient and the corresponding
    index

    Parameters
    ----------
    var : array
        The variable to investigate.
    dx : array
        The grid spacing in x for the different evaluation points
    MXG : [None|int]
        If this is not None, the routine assumes that the variable contains
        ghost points, and MXG will be subtracted from maxGradInd

    Returns
    -------
    maxGrad : float
        The value of the maximum gradient
    maxGradInd : int
        The index at which the maximum gradient can be found at
    """

    ddxVar = DDX(var, dx)

    # nanmax excludes the NaNs
    maxGradInds = np.where(np.abs(ddxVar) == np.nanmax(np.abs(ddxVar)))

    maxGrad = ddxVar[maxGradInds]

    # Take the firs occurence of the x axis
    maxGradInd = int(maxGradInds[1][0])
    if MXG is not None:
        # Subtract MXG as xguards are collected
        maxGradInd -= MXG

    return maxGrad, maxGradInd
#}}}
