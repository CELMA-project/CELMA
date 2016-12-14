#!/usr/bin/env python

"""
Contains function which collects a variable over several output timesteps
"""

import numpy as np
from boutdata import collect

#{{{safeCollect
def safeCollect(*args, **kwargs):
    #{{{docstring
    """
    Wrapper around collect which sets the data immutable

    Parameters
    ----------
    *args : positional arguments
        Positional arguments to collect
    *kwargs : keyword arguments
        Keyword arguments to collect

    Return
    ------
    data : array
        The array is not writeable
    """
    #}}}

    data = collect(*args, **kwargs)
    # write = False prevents writing
    data.setflags(write=False)

    return data
#}}}

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
                    safeCollect(var,\
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

#{{{collectPointTime
def collectPointTime(paths, varName, xInd, yInd, zInd, tInd = None):
    #{{{docstring
    """
    Collects the variable in one spatial point

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : string
        Name of the variable.
    xInd : int
        xInd to collect from
    yInd : int
        yInd to collect from
    zInd : int
        zInd to collect from
    tInd : [None|tuple]
        Start and end of the time if not None

    Returns
    -------
    var : 4d-array
        The time trace of the variable in the position on with shape
        (nt,1,1,1)
    time : 1d-array
        Array of the time at the difference time indices
    """
    #}}}

    if varName == "n":
        collectVarName = "lnN"
    else:
        collectVarName = varName

        varTimeDict = collectiveCollect(paths, [collectVarName, "t_array"],\
                                        collectGhost = False,\
                                        xInd = (xInd, xInd),\
                                        yInd = (yInd, yInd),\
                                        zInd = (zInd, zInd),\
                                        tInd = tInd
                                        )

    var  = varTimeDict[collectVarName]
    time = varTimeDict["t_array"]

    if varName == "n":
        var = np.exp(var)

    return var, time
#}}}

#{{{collectRadialProfileTime
def collectRadialProfileTime(paths, varName, yInd, zInd, tInd = None):
    #{{{docstring
    """
    Collects the variable in one spatial point

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : string
        Name of the variable.
    yInd : int
        yInd to collect from
    zInd : int
        zInd to collect from
    tInd : [None|tuple]
        Start and end of the time if not None

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,nx,1,1)
    time : 1d-array
        Array of the time at the difference time indices
    """
    #}}}

    if varName == "n":
        collectVarName = "lnN"
    else:
        collectVarName = varName

    # Collect the variable
    varTimeDict = collectiveCollect(paths, [collectVarName, "t_array"],\
                                    collectGhost = False,\
                                    yInd = (yInd, yInd),\
                                    zInd = (zInd, zInd),\
                                    tInd = tInd        ,\
                                    )

    var  = varTimeDict[collectVarName]
    time = varTimeDict["t_array"]

    if varName == "n":
        var = np.exp(var)

    return var, time
#}}}

#{{{collectPoloidalProfileTime
def collectPoloidalProfileTime(paths, varName, xInd, yInd, tInd=None):
    #{{{docstring
    """
    Collects the variable in one spatial point

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : string
        Name of the variable.
    xInd : int
        xInd to collect from
    yInd : int
        yInd to collect from
    tInd : [None|tuple]
        Start and end of the time if not None

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,1,1,nz)
    time : 1d-array
        Array of the time at the difference time indices
    """
    #}}}

    if varName == "n":
        collectVarName = "lnN"
    else:
        collectVarName = varName

    # Collect the variable
    varTimeDict = collectiveCollect(paths, [collectVarName, "t_array"],\
                                    collectGhost = False,\
                                    xInd = (xInd, xInd),\
                                    yInd = (yInd, yInd),\
                                    tInd = tInd        ,\
                                    )

    var  = varTimeDict[collectVarName]
    time = varTimeDict["t_array"]

    if varName == "n":
        var = np.exp(var)

    return var, time
#}}}
