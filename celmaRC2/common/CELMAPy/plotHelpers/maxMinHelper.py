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
    vMax : tuple
        The maximum. One per frame
    vMin : tuple
        The minimum. One per frame
    """
    #}}}

    # Get the number of frames
    nFrames = tupleOfArrays[0].shape[0]

    if not(fluct) and not(varyMaxMin):
        # Find the global max and min
        flatArrays = tuple(array.flatten() for array in tupleOfArrays)
        tupleOfArrays = np.concatenate(flatArrays)
        vMax = (np.max(tupleOfArrays),)*nFrames
        vMin = (np.min(tupleOfArrays),)*nFrames
    elif not(fluct) and varyMaxMin:
        # Find the max and min for each frame
        lenTuple = len(tupleOfArrays)
        vMax = np.zeros(nFrames)
        vMin = np.zeros(nFrames)

        for frame in range(nFrames):
            curFrameArrayMax = np.zeros(lenTuple)
            curFrameArrayMin = np.zeros(lenTuple)
            for arrayNr in range(lenTuple):
                curFrameArrayMax[arrayNr] =\
                    np.max(tupleOfArrays[arrayNr][frame,:,:])
                curFrameArrayMin[arrayNr] =\
                    np.min(tupleOfArrays[arrayNr][frame,:,:])

            vMax[frame] = np.max(curFrameArrayMax)
            vMin[frame] = np.min(curFrameArrayMin)

        vMax = tuple(vMax)
        vMin = tuple(vMin)

    elif fluct and not(varyMaxMin):
        # Max and min will be set symmetric
        flatArrays = tuple(array.flatten() for array in tupleOfArrays)
        tupleOfArrays = np.concatenate(flatArrays)

        curMax = np.max(tupleOfArrays)
        curMin = np.min(tupleOfArrays)

        absMax = np.max(np.abs((curMax, curMin)))

        vMax = ( absMax,)*nFrames
        vMin = (-absMax,)*nFrames

    elif fluct and varyMaxMin:
        # Find the max and min for each frame
        lenTuple = len(tupleOfArrays)
        vMax = np.zeros(nFrames)
        vMin = np.zeros(nFrames)

        for frame in range(nFrames):
            curFrameArrayMax = np.zeros(lenTuple)
            curFrameArrayMin = np.zeros(lenTuple)
            for arrayNr in range(lenTuple):
                curFrameArrayMax[arrayNr] =\
                    np.max(tupleOfArrays[arrayNr][frame,:,:])
                curFrameArrayMin[arrayNr] =\
                    np.min(tupleOfArrays[arrayNr][frame,:,:])

            curMax = np.max(curFrameArrayMax)
            curMin = np.min(curFrameArrayMin)

            absMax = np.max(np.abs((curMax, curMin)))

            vMax[frame] = absMax
            vMin[frame] = -absMax

        vMax = tuple(vMax)
        vMin = tuple(vMin)

    return vMax, vMin
#}}}

#{{{getLevelsAnimation
def getLevelsAnimation(vMax, vMin, nCont):
    #{{{docstring
    """
    Sets the contour levels from the tuple of max and min for an animation.

    Parameters
    ----------
    vMax : tuple of nd-arrays
        Tuple of the max. One for each frame.
    vMin : tuple of nd-arrays
        Tuple of the min. One for each frame.
    nCont : int
        Number of contours to use in the contour plot.

    Returns
    -------
    levels : tuple
        The levels to use in the contour plot. One per frame
    """
    #}}}

    nFrames = len(vMax)
    levels = []

    for frame in range(nFrames):
        level = np.linspace(vMin[frame]  ,\
                            vMax[frame]  ,\
                            nCont        ,\
                            endpoint = True)

        if np.min(np.diff(level)) < 0.0:
            levels.append(None)
        else:
            levels.append(level)

    return tuple(levels)
#}}}

#{{{getVmaxVminLevels
def getVmaxVminLevels(plotSuperKwargs,\
                      tupleOfArrays  ,\
                      fluct          ,\
                      varyMaxMin     ,\
                      nLevels = 100  ,\
                      ):
    #{{{docstring
    """
    Getter for vmax, vmin, levels

    Parameters
    ----------
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tupleOfArrays : tuple
        Tuple of the arrays to find the max and min of..
    fluct : bool
        If the variables are fluctuating or not (will symmetrize the max
        and the min).
        Only active if plotSuperKwargs["vmax"] is not set.
    varyMaxMin : bool
        If the max and min allowed to vary
        Only active if plotSuperKwargs["vmax"] is not set.
    nLevels : int
        Number of levels to use.

    Returns
    -------
    vmax : tuple
        Tuple of the vmax values used in contourf.
    vmin : tuple
        Tuple of the vmin values used in contourf.
    levels : tuple
        Tuple of the levels values used in contourf.
    """
    #}}}
    if "vmax" in plotSuperKwargs.keys():
        if plotSuperKwargs["vmax"] is not None:
            vmax = plotSuperKwargs.pop("vmax")
            vmin = plotSuperKwargs.pop("vmin")
        else:
            plotSuperKwargs.pop("vmax")
            plotSuperKwargs.pop("vmin")
            vmax, vmin = getMaxMinAnimation(tupleOfArrays,\
                                            fluct        ,\
                                            varyMaxMin)
    else:
        vmax, vmin = getMaxMinAnimation(tupleOfArrays,\
                                        fluct        ,\
                                        varyMaxMin)
    if "levels" in plotSuperKwargs.keys():
        if plotSuperKwargs["levels"] is not None:
            levels = plotSuperKwargs.pop("levels")
        else:
            plotSuperKwargs.pop("levels")
            levels = getLevelsAnimation(vmax, vmin, nLevels)
    else:
        levels = getLevelsAnimation(vmax, vmin, nLevels)

    return vmax, vmin, levels
#}}}
