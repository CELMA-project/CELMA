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
                    np.gradient(var[tInd, :, yInd, zInd], dx, edge_order=2)

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
                    np.gradient(var[tInd, xInd, :, zInd], dy, edge_order=2)

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
    if len(J.shape) != 2:
       raise ValueError("Input Jacobian must be 2-dimensional")

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
    #{{{docstring
    """
    Returns the index of the maximum absolute gradient along rho.

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
    maxGradInd : int
        The index at which the maximum gradient can be found at
    """
    #}}}

    ddxVar = DDX(var, dx)

    # nanmax excludes the NaNs
    maxGradInds = np.where(np.abs(ddxVar) == np.nanmax(np.abs(ddxVar)))

    maxGrad = ddxVar[maxGradInds]

    # Take the firs occurence of the x axis
    maxGradInd = int(maxGradInds[1][0])
    if MXG is not None:
        # Subtract MXG as xguards are collected
        maxGradInd -= MXG

    return maxGradInd
#}}}

#{{{findLargestParallelGrad
def findLargestParallelGrad(var, dy, MYG = None):
    #{{{docstring
    """
    Returns the index of the maximum absolute gradient along z.

    Parameters
    ----------
    var : array
        The variable to investigate.
    dy : array
        The grid spacing in y for the different evaluation points
    MYG : [None|int]
        If this is not None, the routine assumes that the variable contains
        ghost points, and MYG will be subtracted from maxGradInd

    Returns
    -------
    maxGradInd : int
        The index at which the maximum gradient can be found at
    """
    #}}}

    ddyVar = DDY(var, dy)

    # nanmax excludes the NaNs
    maxGradInds = np.where(np.abs(ddyVar) == np.nanmax(np.abs(ddyVar)))

    maxGrad = ddyVar[maxGradInds]

    # Take the firs occurence of the y axis
    maxGradInd = int(maxGradInds[2][0])
    if MYG is not None:
        # Subtract MYG as yguards are collected
        maxGradInd -= MYG

    return maxGradInd
#}}}

#{{{findLargestPoloidalGrad
def findLargestPoloidalGrad(var, J):
    #{{{docstring
    """
    Returns the index of the maximum absolute gradient along theta

    Parameters
    ----------
    var : array
        The variable to investigate.
    J : array
        The Jacobian

    Returns
    -------
    maxGradInd : int
        The index at which the maximum gradient can be found at
    """
    #}}}

    ddzVar = DDZ(var, J)

    # nanmax excludes the NaNs
    maxGradInds = np.where(np.abs(ddzVar) == np.nanmax(np.abs(ddzVar)))

    maxGrad = ddzVar[maxGradInds]

    # Take the firs occurence of the z axis
    maxGradInd = int(maxGradInds[3][0])

    return maxGradInd
#}}}

#{{{findLargestRadialGradN
def findLargestRadialGradN(steadyStatePath):
    #{{{docstring
    """
    Find the largest gradient in n.

    Parameters
    ----------
    steadyStatePath : str
        Path to collect from

    Returns
    -------
    xInd : int
        Index of largest n
    """
    #}}}

    dx = getUniformSpacing(steadyStatePath, "x")
    # Check last t index
    with DataFile(os.path.join(steadyStatePath, "BOUT.dmp.0.nc")) as f:
        tLast = len(f.read("t_array")) - 1

    # In the steady state, the max gradient in "n" is the same
    # throughout in the domain, so we use yInd=0, zInd=0 in the
    # collect
    lnN = collect("lnN",\
                  path=steadyStatePath,\
                  xguards=False,\
                  yguards=False,\
                  tind   = [tLast, tLast],\
                  info=False)
    n = calcN(lnN)
    xInd  = findLargestRadialGrad(n, dx[0,0])

    return xInd
#}}}
