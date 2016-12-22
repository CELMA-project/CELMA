#!/usr/bin/env python

"""
Contains function which collects a variable over several output timesteps
"""

from boutdata import collect
from boututils.datafile import DataFile
import numpy as np
import os

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
                      zInd         = None ):
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
    tInd : [None|tuple]
        Start and end of the time if not None
    xInd : [None|2d array]
        x index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)
    yInd : [None|2d array]
        y index range to collect. The first index is the start, and the
        second is the end of the range (inclusive)
    zInd : [None|2d array]
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

    if tInd is not None and len(paths) > 1:
        # Check if it necessary to use all paths
        paths = removePathsOutsideRange(paths, tInd)
        if len(paths) != 1:
            postCollectTInd = tInd
            tInd = None
        else:
            postCollectTInd = None

    for path in paths:
        for var in varStrings:
            try:
                # Make a local var which is reused for every interation,
                # then concatenate the dictionary
                curVar =\
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

            # Ensure 4D
            if len(curVar.shape) == 3:
                # Make it a 4d variable
                time = collectTime(paths, tInd)
                tmp = np.zeros((len(time), *curVar.shape))
                # Copy the field in to each time
                tmp[:] = curVar
                curVar = tmp

            if tInd is not None and len(paths) > 1:
                # Do slicing post collection if necessary
                if postCollectTInd is not None:
                    curVar = curVar[tInd[0]:tInd[1]]

            # Set data[var] to curVar the first time in order to get the
            # correct dimensions
            if data[var] is None:
                data[var] = curVar
            else:
                # Remove first point in time in the current time as this
                # is the same as the last of the previous
                data[var]=np.concatenate((data[var], curVar[1:,:,:,:]), axis=0)

    return data
#}}}

#{{{removePathsOutsideRange
def removePathsOutsideRange(paths, tInd):
    #{{{docstring
    """
    Removes paths which are outside of the time range.

    Parameters
    ----------
    paths : tuple
        The paths to collect from
    tInd : tuple
        The time indices


    Returns
    -------
    paths : tuple
        The filtered paths
    """
    #}}}
    # Check if it necessary to use all paths
    lenT = 0
    # Cast to list
    paths = list(paths)
    removePaths = []
    # Check if the paths are needed from the rigth
    for ind in range(len(paths)):
        with DataFile(os.path.join(paths[ind],"BOUT.dmp.0.nc")) as f:
            t = f.read("t_array")
            lenT += len(t)
            if tInd[1] < lenT:
                lastInd = ind
            else:
                break
    # Remove from the rigth
    paths = paths[:lastInd+1]
    # Check if the paths are needed from the left
    lenT = 0
    for path in paths:
        with DataFile(os.path.join(path,"BOUT.dmp.0.nc")) as f:
            t = f.read("t_array")
            lenT += len(t)
            if tInd[0] > lenT:
                removePaths.append(removePaths)
            else:
                break

    for remove in removePaths:
        paths.remove(remove)

    # Cast to tuple
    paths = tuple(path)

    return paths
#}}}

#{{{collectTime
def collectTime(paths, tInd = None):
    #{{{docstring
    """
    Collects the time

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    tInd : [None|tuple]
        Start and end of the time if not None

    Returns
    -------
    time : 1d-array
        Array of the time at the difference time indices
    """
    #}}}

    # Initialize
    time = None

    for path in paths:
        with DataFile(os.path.join(path,"BOUT.dmp.0.nc")) as f:
            if time is None:
                time = f.read("t_array")
            else:
                # Remove first point in time in the current time as this
                # is the same as the last of the previous
                time = np.concatenate((time, f.read("t_array")[1:]), axis=0)

    if tInd is not None:
        # +1 as the last point is excluded from a slice
        time = time[tInd[0]:tInd[1]+1]

    return time
#}}}

#{{{collectPoint
def collectPoint(paths, varName, xInd, yInd, zInd, tInd = None):
    #{{{docstring
    """
    Collects the variable in one spatial point

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
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
    """
    #}}}

    # Cast to tuple
    varDict = collectiveCollect(paths, (varName,)   ,\
                                collectGhost = False,\
                                xInd = (xInd, xInd) ,\
                                yInd = (yInd, yInd) ,\
                                zInd = (zInd, zInd) ,\
                                tInd = tInd
                               )

    var  = varDict[varName]

    return var
#}}}

#{{{collectRadialProfile
def collectRadialProfile(paths, varName, yInd, zInd,\
                         tInd = None, collectGhost = False):
    #{{{docstring
    """
    Collects the variable in along a radial line

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
    yInd : tuple
        yInd start and yInd end to collect from
    zInd : tuple
        zInd start and zInd end to collect from
    tInd : [None|tuple]
        Start and end of the time if not None
    collectGhost : bool
        Whether or not the ghost should be collected

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,nx,1,1)
    """
    #}}}

    # Collect the variable
    varDict = collectiveCollect(paths, (varName,)           ,\
                                 collectGhost = collectGhost,\
                                 yInd = yInd                ,\
                                 zInd = zInd                ,\
                                 tInd = tInd                ,\
                                )

    var  = varDict[varName]

    return var
#}}}

#{{{collectParallelProfile
def collectParallelProfile(paths, varName, xInd, zInd,\
                           tInd = None, collectGhost = False):
    #{{{docstring
    """
    Collects the variable in along a parallel line

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
    xInd : tuple
        xInd start and xInd end to collect from
    zInd : tuple
        zInd start and zInd end to collect from
    tInd : [None|tuple]
        Start and end of the time if not None
    collectGhost : bool
        Whether or not the ghost should be collected

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,1,ny,1)
    """
    #}}}

    # Collect the variable
    varDict = collectiveCollect(paths, (varName,)          ,\
                                collectGhost = collectGhost,\
                                xInd = xInd                ,\
                                zInd = zInd                ,\
                                tInd = tInd                ,\
                               )

    var  = varDict[varName]

    return var
#}}}

#{{{collectPoloidalProfile
def collectPoloidalProfile(paths, varName, xInd, yInd,\
                           tInd = None, collectGhost = False):
    #{{{docstring
    """
    Collects the variable along a poloidal line

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
    xInd : tuple
        xInd start and xInd end to collect from
    yInd : tuple
        yInd start and yInd end to collect from
    tInd : [None|tuple]
        Start and end of the time if not None
    collectGhost : bool
        Whether or not the ghost should be collected

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,1,1,nz)
    """
    #}}}

    # Collect the variable
    varDict = collectiveCollect(paths, (varName,)          ,\
                                collectGhost = collectGhost,\
                                xInd = xInd                ,\
                                yInd = yInd                ,\
                                tInd = tInd                ,\
                               )

    var  = varDict[varName]

    return var
#}}}

#{{{collectConstRho
def collectConstRho(paths, varName, xInd,\
                    tInd = None, collectGhost = False):
    #{{{docstring
    """
    Collects the variable with at specified rho

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
    xInd : tuple
        xInd start and xInd end to collect from
    tInd : [None|tuple]
        Start and end of the time if not None
    collectGhost : bool
        Whether or not the ghost should be collected

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,1,ny,nz)
    """
    #}}}

    # Collect the variable
    varDict = collectiveCollect(paths, (varName,)          ,\
                                collectGhost = collectGhost,\
                                xInd = xInd                ,\
                                tInd = tInd                ,\
                               )

    var  = varDict[varName]

    return var
#}}}

#{{{collectConstZ
def collectConstZ(paths, varName, yInd,\
                  tInd = None, collectGhost = False):
    #{{{docstring
    """
    Collects the variable with at specified z

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
    yInd : tuple
        yInd start and yInd end to collect from
    tInd : [None|tuple]
        Start and end of the time if not None
    collectGhost : bool
        Whether or not the ghost should be collected

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,nx,1,nz)
    """
    #}}}

    # Collect the variable
    varDict = collectiveCollect(paths, (varName,)          ,\
                                collectGhost = collectGhost,\
                                yInd = yInd                ,\
                                tInd = tInd                ,\
                               )

    var  = varDict[varName]

    return var
#}}}

#{{{collectConstTheta
def collectConstTheta(paths, varName, zInd,\
                      tInd = None, collectGhost = False):
    #{{{docstring
    """
    Collects the variable with at a specified theta

    Parameters
    -----------
    paths : iterable of strings
        What path to use when collecting the variable. Must be in
        ascending temporal order as the variable will be
        concatenated.
    varName : str
        Name of the variable to collect.
    zInd : tuple
        zInd start and zInd end to collect from
    tInd : [None|tuple]
        Start and end of the time if not None
    collectGhost : bool
        Whether or not the ghost should be collected

    Returns
    -------
    var : 4d-array
        The variable with the shape (nt,nx,ny,1)
    """
    #}}}

    # Collect the variable
    varDict = collectiveCollect(paths, (varName,)          ,\
                                collectGhost = collectGhost,\
                                zInd = zInd                ,\
                                tInd = tInd                ,\
                               )

    var  = varDict[varName]

    return var
#}}}
