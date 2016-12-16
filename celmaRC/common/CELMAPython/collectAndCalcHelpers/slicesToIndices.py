#!/usr/bin/env python

"""
Contains function which converts a slice to indices used for BOUT++ collection.
"""

import numpy as np

#{{{slicesToIndices
def slicesToIndices(path, theSlice, dimension, xguards=False, yguards=False):
    #{{{docstring
    """
    Return the slice such that it can be given as an input to "collect"

    Parameters
    ----------
    theSlice : [slice | int | None]
        Current slice to use
    dimension : ["x" | "y" | "z" | "t"]
        The dimension to slice in

    Returns
    -------
    indices : tuple
        Tuple containing the start and the stop values from the slice
    """
    #}}}

    if type(theSlice) == slice:
        indices = []
        indices.append(theSlice.start)
        if theSlice.stop == None:
            # Find the last index
            if dimension == "x":
                dimLen = getSizes(path, dimension, includeGhost = xguards)
            elif dimension == "y":
                dimLen = getSizes(path, dimension, includeGhost = yguards)
            elif dimension == "z":
                dimLen = getSizes(path, dimension)
            elif dimension == "t":
                dimLen = getSizes(path, dimension)
            else:
                raise ValueError("Unknown coordinate {}".format(coordinate))

            # Subtract 1 in the end as indices counts from 0
            indices.append(dimLen - 1)
            indices = tuple(indices)
        else:
            indices.append(theSlice.stop)

    elif type(theSlice) is int:
        indices = (theSlice, theSlice)
    elif theSlice is None:
        indices = theSlice
    else:
        message = "{} not allowed as an input slice".format(type(theSlice))
        raise ValueError(message)

    # Check for negative indices
    if indices is not None:
        for ind in range(len(indices)):
            if indices[ind] < 0:
                # Find the last index
                if dimension == "x":
                    dimLen = getSizes(path, dimension, includeGhost = xguards)
                elif dimension == "y":
                    dimLen = getSizes(path, dimension, includeGhost = yguards)
                elif dimension == "z":
                    dimLen = getSizes(path, dimension)
                elif dimension == "t":
                    dimLen = getSizes(path, dimension)
                else:
                    raise ValueError("Unknown coordinate {}".format(coordinate))

                # Subtract 1 in the end as indices counts from 0
                realInd = dimLen + indices[ind] - 1

                if realInd < 0:
                    message  = ("Index {0} out of range for {1}"
                                ", as {1} has only {2} elements").\
                        format(indices[ind], dimension, dimLen)
                    raise IndexError(message)
            if ind == 0:
                start = realInd
            else:
                end = realInd

        # Cast to tuple
        indices = tuple(start, end)

    return indices
#}}}
