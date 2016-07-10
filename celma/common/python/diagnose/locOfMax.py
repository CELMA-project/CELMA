#!/usr/bin/env python

import numpy as np

def locOfMax(f):
    """
    Finds the location of the max abs occurence of a field

    Parameters
    ----------
    f : BOUT++ field
        The field to find the location from

    Returns
    loc: array
        Array of the location of the max abs occurence
    -------
    """

    fAMax = np.abs(f.max())
    fMin  = f.min()
    fAMin = np.abs(f.min())
    if fAMax > fAMin:
        print("Max abs value is positive")
        loc = np.where(np.isclose(f, fAMax))
    else:
        print("Max abs value is negative")
        loc = np.where(np.isclose(f, fMin))

    return loc
