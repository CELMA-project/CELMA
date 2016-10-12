#!/usr/bin/env python

"""
Contains function which collects a variable over several output timesteps
"""

import numpy as np
from boutdata import collect

#{{{collectiveCollect
def collectiveCollect(paths               ,\
                      varStrings          ,\
                      collectGhost = False,\
                      tInd         = None ,\
                      yInd         = None ,\
                      xInd         = None ,\
                      zInd         = None):
    #{{{docstring
    """
    Collects variables from several paths

    Parameters
    ----------
    paths : iterable of strings
        The paths to collect from. Must be in ascending order of the
        simulation time, as the variables are being concatenated
    varStrings : iterable of strings
        The variables to be collected
    collectGhost : bool
        If the ghost is to be collected
    tind : [None|2d array]
        Time index range to collect. The first index is the start, and
        the second is the end of the range (inclusive)
    xind : [None|2d array]
        x index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)
    yind : [None|2d array]
        y index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)
    zind : [None|2d array]
        z index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)

    Return
    ------
    data : dict
        A dictionary of the concatenated variables
    """
    #}}}

    # Initialize the data
    data = {var: None for var in varStrings}

    for path in paths:
        for var in varStrings:

            try:
                # Make a local var which is reused for every interation,
                # then concatenate the dictionary
                localVar =\
                    collect(var,\
                            path     = path        ,\
                            tind     = tInd        ,\
                            xind     = xInd        ,\
                            yind     = yInd        ,\
                            zind     = zInd        ,\
                            xguards  = collectGhost,\
                            yguards  = collectGhost,\
                            info     = False        )
            except OSError:
                # An OSError is thrown if the file is not found
                raise ValueError("No collectable files found in {}".\
                                 format(path))


            # Set data[var] to localVar the first time in order to get the
            # correct dimensions
            if data[var] is None:
                data[var] = localVar
            else:
                data[var] = np.concatenate((data[var], localVar), axis=0)

    return data
#}}}
