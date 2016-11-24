#!/usr/bin/env python

"""
Contains functions for taking poloidal averages
"""

import numpy as np

#{{{polAvg
def polAvg(f):
    """
    Returns the poloidal average of a field.

    Parameters
    ----------
    f : array
        The field to find the poloidal average of.
        The field must be a 4D field, and should not include the last
        poloidal slice (i.e. the domain should go from [0,2pi[)

    Returns
    -------
    out : array
        The poloidal average of the field
    """

    tLen, xLen, yLen, zLen = f.shape
    out = np.zeros(f.shape)

    for t in range(tLen):
        for x in range(xLen):
            for y in range(yLen):
                out[t,x,y,:] = f[t,x,y,:].mean()

    return out
#}}}

#{{{timeAvg
def timeAvg(f, t, startInd = 0, endInd = -1):
    """
    Returns the poloidal average of a field.

    Parameters
    ----------
    f : array
        The field to find the time average of.
        The field must be a 4D field.
    t : array
        The time.
        Must have the same temporal dimension as f.
    startInd : int
        Start index to take the average from
    endInd : int
        End index to take the average from
    stdDev : bool
        Whether or not to calculate the standard deviation and include
        this in the output.

    Returns
    -------
    NOTE: Output arrays will have time dimensions of
    np.floor(endInd - startInd)/len(t)

    out : array
        The poloidal average of the field
    t : array
        Averaged time array
    """

    tLen, xLen, yLen, zLen = f.shape

    # Check for negative indices
    if startInd < 0:
        # Subtract 1 as indices count from 0
        startInd = (tLen - 1) + startInd
        if startInd < 0:
            raise ValueError("startInd={} is out of range".\
                    format(startInd))
    if endInd < 0:
        # Subtract 1 as indices count from 0
        endInd = (tLen - 1) + endInd
        if endInd < 0:
            raise ValueError("endInd={} is out of range".\
                    format(endInd))

    diffT = endInd - startInd
    if dffT < 0:
        raise ValueError("Must have endInd>startInd")

    tLenOut = np.floor(endInd - startInd)/tLen
    outDim  = (tLenOut, f.shape[1], f.shape[2], f.shape[3])

    out = np.zeros(outDim)

    for avgTInd in range(tLenOut):
        tStart = startInd + avgTInd
        tEnd   = endInd   + avgTInd
        for x in range(xLen):
            for y in range(yLen):
                for z in range(zLen):
                    out[avgTInd,x,y,z] = f[tStart:tEnd,x,y,z].mean()

    return out, t[startInd:startInd + tLenOut]
#}}}
