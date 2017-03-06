#!/usr/bin/env python

"""
Contains functions dealing with sizes of the time
"""

from boututils.datafile import DataFile
import os

#{{{getTSize
def getTSize(paths):
    #{{{docstring
    """
    Fastest way to obtain the time size.

    Parameters
    ----------
    paths : tuple
        Tuple of paths to collect from

    Returns
    -------
    tSize : int
        Size of the time
    """
    #}}}

    tSize = 0
    for path in paths:
        with DataFile(os.path.join(path,"BOUT.dmp.0.nc")) as f:
            t = f.read("t_array")
            tSize += len(t)

    return tSize
#}}}
