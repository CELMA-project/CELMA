#!/usr/bin/env python

"""
Contains functions for taking poloidal averages
"""

import numpy as np

#{{{polAvg
def polAvg(f):
    #{{{docstring
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
    #}}}

    tLen, xLen, yLen, zLen = f.shape
    out = np.zeros(f.shape)

    for t in range(tLen):
        for x in range(xLen):
            for y in range(yLen):
                out[t,x,y,:] = f[t,x,y,:].mean()

    return out
#}}}

#{{{timeAvg
def timeAvg(f, t = None, startInd = 0, endInd = -1):
    #{{{docstring
    """
    Returns the poloidal average of a field.

    The average is a sliding the time-window through the time such that
    The output arrays will have time dimensions of
    int(np.floor((tLen-1)/(endInd - startInd)))

    Parameters
    ----------
    f : array
        The field to find the time average of.
        The field must be a 4D field.
    t : [None|array]
        The time.
        Must have the same temporal dimension as f.
    startInd : int
        Start index to take the average from.
    endInd : int
        End index to take the average from.

    Returns
    -------
    avgF : array
        The poloidal average of the field.
    avgT : array
        Only an output if t is not None.
        Averaged time array.
    """
    #}}}

    tLen, xLen, yLen, zLen = f.shape

    # Check for negative indices
    if startInd < 0:
        # Subtract 1 as indices count from 0
        startInd = (tLen - 1) + startInd
        if startInd < 0:
            raise ValueError("startInd={} is out of range".\
                    format(startInd))
    if endInd < 0:
        endInd = tLen + endInd
        if endInd < 0:
            raise ValueError("endInd={} is out of range".\
                    format(endInd))

    diffT = endInd - startInd
    if diffT < 0:
        raise ValueError("Must have endInd>startInd")

    tLenOut = int(np.floor((tLen-1)/(endInd - startInd)))
    outDim  = (tLenOut, xLen, yLen, zLen)

    outF = np.zeros(outDim)

    if t is not None:
        outT = np.zeros(tLenOut)

        # As the average will be sliding, we now the averaged times in advance
        startTAvgInd = int(t[startInd:endInd+1].mean())
        endTAvgInd   = startTAvgInd + tLenOut

        outT = np.array(range(startTAvgInd, endTAvgInd))

    for avgTInd in range(tLenOut):
        tStart = startInd + avgTInd
        tEnd   = endInd   + avgTInd
        for x in range(xLen):
            for y in range(yLen):
                for z in range(zLen):
                    # +1 as slicing does not include last point
                    outF[avgTInd,x,y,z] = f[tStart:tEnd+1,x,y,z].mean()

    if t is not None:
        return outF, outT
    else:
        return outF
#}}}
