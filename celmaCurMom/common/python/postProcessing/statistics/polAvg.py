#!/usr/bin/env python

"""
Contains function for taking poloidal average
"""

#{{{polAvg
def polAvg(f):
    """
    Takes the average poloidally.

    Input
    f   - The field before taking the poloidal average
          The field must be a 4D field, and should not include the last
          poloidal slice.

    Output
    f   - The field after taking the poloidal average
    """

    tLen, xLen, yLen, zLen = f.shape

    for t in range(tLen):
        for x in range(xLen):
            for y in range(yLen):
                f[t,x,y,:] -= f[t,x,y,:].mean()

    return f
#}}}
