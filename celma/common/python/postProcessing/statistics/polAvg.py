#!/usr/bin/env python

"""
Contains function for taking poloidal average
"""

import numpy as np

#{{{polAvg
def polAvg(f):
    """
    Takes the average poloidally.

    Input
    f   - The field before taking the poloidal average
          The field must be a 4D field, and should not include the last
          poloidal slice.

    Output
    out   - The field after taking the poloidal average
    """

    tLen, xLen, yLen, zLen = f.shape
    out = np.zeros(f.shape)

    for t in range(tLen):
        for x in range(xLen):
            for y in range(yLen):
                out[t,x,y,:] = f[t,x,y,:] - f[t,x,y,:].mean()

    return out
#}}}
