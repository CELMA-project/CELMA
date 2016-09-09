#!/usr/bin/env python

"""
Contains function for taking poloidal average
"""

import numpy as np

#{{{polAvg
def polAvg(f):
    """
    Returns the poloidal average of a field.

    Input
    f   - The field to find the poloidal average of.
          The field must be a 4D field, and should not include the last
          poloidal slice (i.e. the domain should go from [0,2pi[)

    Output
    out   - The poloidal average of the field
    """

    tLen, xLen, yLen, zLen = f.shape
    out = np.zeros(f.shape)

    for t in range(tLen):
        for x in range(xLen):
            for y in range(yLen):
                out[t,x,y,:] = f[t,x,y,:].mean()

    return out
#}}}
