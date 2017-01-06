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

    # Initialize postCollectTInd for later checks
    postCollectTInd = None
    if len(paths) > 1:
        if tInd is not None:
            # Check if it necessary to use all paths
            paths, tInd = removePathsOutsideRange(paths, tInd)
            # If more than one path, we will slice in time after
            # collection
            postCollectTInd = tInd
            tInd = None
        else:
            postCollectTInd = None

    # Must cast to list due to corner cases in collect
    if tInd is not None:
        tInd = list(tInd)

    for var in varStrings:
        for path in paths:
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
            elif len(curVar.shape) == 1:
                # Make it a 4d variable
                tmp = np.zeros((len(curVar), 1, 1, 1))
                # Copy the field in to each time
                tmp[:,0,0,0] = curVar
                curVar = tmp

            # Set data[var] to curVar the first time in order to get the
            # correct dimensions
            if data[var] is None:
                data[var] = curVar
            else:
                # Remove first point in time in the current time as this
                # is the same as the last of the previous
                data[var]=np.concatenate((data[var], curVar[1:,:,:,:]), axis=0)

        # Do slicing post collection if necessary
        if postCollectTInd is not None:
            # +1 to include the last point
            lastPost = postCollectTInd[1]+1 if postCollectTInd[1] is not None\
                       else None

            data[var] = data[var][postCollectTInd[0]:lastPost]

    # Recast to tuple
    if tInd is not None:
        tInd = tuple(tInd)

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
    tInd : tuple
        Updated time indices
    """
    #}}}

# FIXME:
    from copy import copy
    startPaths = copy(paths)
    startTInd  = copy(tInd)
# END FIXME
    # Cast to list
    paths = list(paths)

    # Check if the paths are needed from the rigth
    if tInd[1] is not None:
        lenT = 0

        lastInd = None
        for ind in range(len(paths)):
            with DataFile(os.path.join(paths[ind],"BOUT.dmp.0.nc")) as f:
                t = f.read("t_array")
                lenT += len(t)
                if tInd[1] < lenT:
                    lastInd = ind+1
                    break

        # Remove from the rigth
        paths = paths[:lastInd]

    # Check if the paths are needed from the left
    if tInd[0] is not None:
        lenT  = 0
        # If folders are removed, subtract tInd with number of removed indices
        # lenTs is a placeholder for this
        removePaths = []
        lenTs       = []
        for path in paths:
            with DataFile(os.path.join(path,"BOUT.dmp.0.nc")) as f:
                t = f.read("t_array")
                curLenT = len(t)
                lenTs.append(curLenT)
                lenT += curLenT
                if tInd[0] > lenT:
                    removePaths.append(path)
                else:
                    break

        tInd = list(tInd)
        for remove, lenT in zip(removePaths, lenTs):
            # +1 removes the duplicate which is present in both files
#            tInd[0] = (tInd[0] - lenT) + 1
            tInd[0] = (tInd[0] - lenT)
            tInd[1] =  tInd[1] - lenT if tInd[1] is not None else None
            paths.remove(remove)
        # Remove the duplicates present in more files by shifting
        # the first time index by one per duplicate
        tInd[0] += len(paths)-1
        tInd = tuple(tInd)

    # Cast to tuple
# FIXME:
    inds = [startPaths.index(path) for path in paths if path in startPaths]
    if tInd[0] is None:
        b = "tInd[0] is None"
        if tInd[1] is None:
            b = "tInd is None"
    elif tInd[1] is None:
        b = "tInd[1] is None"
    else:
        print("tInd[0]-tInd[1]={}".format(tInd[0]-tInd[1]))
        b = (tInd[0]-tInd[1]) == (startTInd[0]-startTInd[1])
    print("\nstartPathLen = {}, endPaths = {}, tInd==startTInd = {}\n".\
            format(len(startPaths),  inds, b)
            )
    # raise RuntimeError("YOU QUITTED")
    # import pdb; pdb.set_trace()
# END FIXME:
    paths = tuple(paths)

    return paths, tInd
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
        lastT = tInd[1]+1 if tInd[1] is not None else None
        time  = time[tInd[0]:lastT]

# FIXME
    print("len(time) = {}".format(len(time)))
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
