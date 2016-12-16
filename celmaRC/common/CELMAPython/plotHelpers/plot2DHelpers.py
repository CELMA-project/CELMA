#!/usr/bin/env python

"""
Contains helper functions for 2D plots
"""

#{{{getMaxMinAnimation
def getMaxMinAnimation(tupleOfArrays, nFrames, fluct, varyMaxMin):
    #{{{docstring
    """
    Finds the max and min for each frame in the animation.

    Parameters
    ----------
    tupleOfArrays : tuple of nd-arrays
        Tuple of the arrays to find the max and min of.
    nFrames : int
        Number of frames.
    fluct : bool
        Whether or not the max and min should be symmetric around 0.
    varyMaxMin : bool
        Whether or not the max and min are allowed to vary from frame to
        frame.

    Returns
    -------
    varMax : tuple
        The maximum. One per frame
    varMin : tuple
        The minimum. One per frame
    """
    #}}}

    if not(fluct) and not(varyMaxMin):
        # Find the global max and min
        flatArrays = tuple(array.flatten() for array in tupleOfArrays)
        tupleOfArrays = np.concatenate(flatArrays)
        varMax = (np.max(tupleOfArrays),)*nFrames
        varMin = (np.min(tupleOfArrays),)*nFrames
    elif not(fluct) and varyMaxMin:
        # Find the max and min for each frame
        lenTuple = len(tupleOfArrays)
        varMax = np.zeros(nFrames)
        varMin = np.zeros(nFrames)
        tupleOfMax = np.zeros(lenTuple)
        tupleOfMin = np.zeros(lenTuple)

        for frame in range(nFrames):
            curFrameArrayMax = np.zeros(nFrames)
            curFrameArrayMin = np.zeros(nFrames)
            for arrayNr in range(lenTuple):
                curFrameArrayMax[arrayNr] = np.max(tupleOfArrays[arrayNr])
                curFrameArrayMin[arrayNr] = np.min(tupleOfArrays[arrayNr])

            varMax[frame] = np.max(curFrameArrayMax)
            varMin[frame] = np.min(curFrameArrayMin)

        varMax = tuple(varMax)
        varMin = tuple(varMin)

    elif fluct and not(varyMaxMin):
        # Max and min will be set symmetric
        flatArrays = tuple(array.flatten() for array in tupleOfArrays)
        tupleOfArrays = np.concatenate(flatArrays)

        curMax = np.max(tupleOfArrays)
        curMin = np.min(tupleOfArrays)

        absMax = np.max(np.abs(curMax, curMin))

        varMax = ( absMax,)*nFrames
        varMin = (-absMax,)*nFrames

    elif fluct and varyMaxMin:
        # Find the max and min for each frame
        lenTuple = len(tupleOfArrays)
        varMax = np.zeros(nFrames)
        varMin = np.zeros(nFrames)
        tupleOfMax = np.zeros(lenTuple)
        tupleOfMin = np.zeros(lenTuple)

        for frame in range(nFrames):
            curFrameArrayMax = np.zeros(nFrames)
            curFrameArrayMin = np.zeros(nFrames)
            for arrayNr in range(lenTuple):
                curFrameArrayMax[arrayNr] = np.max(tupleOfArrays[arrayNr])
                curFrameArrayMin[arrayNr] = np.min(tupleOfArrays[arrayNr])

            curMax = np.max(curFrameArrayMax)
            curMin = np.max(curFrameArrayMin)

            absMax = np.max(np.abs(curMax, curMin))

            varMax[frame] = absMax
            varMin[frame] = -absMax

        varMax = tuple(varMax)
        varMin = tuple(varMin)

    return varMax, varMin
#}}}

#{{{getLevelsAnimation
def getLevelsAnimation(varMax, varMin, nCont):
    #{{{docstring
    """
    Sets the contour levels from the tuple of max and min for an animation.

    Parameters
    ----------
    varMax : tuple of nd-arrays
        Tuple of the max. One for each frame.
    varMin : tuple of nd-arrays
        Tuple of the min. One for each frame.
    nCont : int
        Number of contours to use in the contour plot.

    Returns
    -------
    levels : tuple
        The levels to use in the contour plot. One per frame
    """
    #}}}

    nFrames = len(varMax)
    levels = np.zeros(nFrames)

    for frame in range(nFrames):
        level = np.linspace(varMin  ,\
                            varMax  ,\
                            nCont   ,\
                            endpoint = True)

        if np.amin(np.diff(levels)) > 0.0:
            levels[frame] = None
        else:
            levels[frame] = level

    return tuple(levels)
#}}}
