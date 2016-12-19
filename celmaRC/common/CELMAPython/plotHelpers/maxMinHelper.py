#!/usr/bin/env python

"""
Contains helper functions for setting plot limits
"""

import numpy as np

#{{{getMaxMinAnimation
def getMaxMinAnimation(tupleOfArrays, fluct, varyMaxMin):
    #{{{docstring
    """
    Finds the max and min for each frame in the animation.

    Parameters
    ----------
    tupleOfArrays : tuple of nd-arrays
        Tuple of the arrays to find the max and min of.
    fluct : bool
        Whether or not the max and min should be symmetric around 0.
    varyMaxMin : bool
        Whether or not the max and min are allowed to vary from frame to
        frame.

    Returns
    -------
    vmax : tuple
        The maximum. One per frame
    vmin : tuple
        The minimum. One per frame
    """
    #}}}

    # Get the number of frames
    nFrames = tupleOfArrays[0].shape[0]

    if not(fluct) and not(varyMaxMin):
        # Find the global max and min
        flatArrays = tuple(array.flatten() for array in tupleOfArrays)
        tupleOfArrays = np.concatenate(flatArrays)
        vmax = (np.max(tupleOfArrays),)*nFrames
        vmin = (np.min(tupleOfArrays),)*nFrames
    elif not(fluct) and varyMaxMin:
        # Find the max and min for each frame
        lenTuple = len(tupleOfArrays)
        vmax = np.zeros(nFrames)
        vmin = np.zeros(nFrames)

        for frame in range(nFrames):
            curFrameArrayMax = np.zeros(nFrames)
            curFrameArrayMin = np.zeros(nFrames)
            for arrayNr in range(lenTuple):
                curFrameArrayMax[arrayNr] =\
                    np.max(tupleOfArrays[arrayNr][frame,:,:])
                curFrameArrayMin[arrayNr] =\
                    np.min(tupleOfArrays[arrayNr][frame,:,:])

            vmax[frame] = np.max(curFrameArrayMax)
            vmin[frame] = np.min(curFrameArrayMin)

        vmax = tuple(vmax)
        vmin = tuple(vmin)

    elif fluct and not(varyMaxMin):
        # Max and min will be set symmetric
        flatArrays = tuple(array.flatten() for array in tupleOfArrays)
        tupleOfArrays = np.concatenate(flatArrays)

        curMax = np.max(tupleOfArrays)
        curMin = np.min(tupleOfArrays)

        absMax = np.max(np.abs((curMax, curMin)))

        vmax = ( absMax,)*nFrames
        vmin = (-absMax,)*nFrames

    elif fluct and varyMaxMin:
        # Find the max and min for each frame
        lenTuple = len(tupleOfArrays)
        vmax = np.zeros(nFrames)
        vmin = np.zeros(nFrames)

        for frame in range(nFrames):
            curFrameArrayMax = np.zeros(nFrames)
            curFrameArrayMin = np.zeros(nFrames)
            for arrayNr in range(lenTuple):
                curFrameArrayMax[arrayNr] =\
                    np.max(tupleOfArrays[arrayNr][frame,:,:])
                curFrameArrayMin[arrayNr] =\
                    np.min(tupleOfArrays[arrayNr][frame,:,:])

            curMax = np.max(curFrameArrayMax)
            curMin = np.min(curFrameArrayMin)

            absMax = np.max(np.abs((curMax, curMin)))

            vmax[frame] = absMax
            vmin[frame] = -absMax

        vmax = tuple(vmax)
        vmin = tuple(vmin)

    return vmax, vmin
#}}}

#{{{getLevelsAnimation
def getLevelsAnimation(vmax, vmin, nCont):
    #{{{docstring
    """
    Sets the contour levels from the tuple of max and min for an animation.

    Parameters
    ----------
    vmax : tuple of nd-arrays
        Tuple of the max. One for each frame.
    vmin : tuple of nd-arrays
        Tuple of the min. One for each frame.
    nCont : int
        Number of contours to use in the contour plot.

    Returns
    -------
    levels : tuple
        The levels to use in the contour plot. One per frame
    """
    #}}}

    nFrames = len(vmax)
    levels = []

    for frame in range(nFrames):
        level = np.linspace(vmin[frame]  ,\
                            vmax[frame]  ,\
                            nCont        ,\
                            endpoint = True)

        if np.min(np.diff(level)) < 0.0:
            levels.append(None)
        else:
            levels.append(level)

    return tuple(levels)
#}}}
