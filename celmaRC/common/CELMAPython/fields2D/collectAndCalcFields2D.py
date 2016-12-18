#!/usr/bin/env python

"""
Contains function for collecting and calculating the 2D fields
"""

from ..unitsConverter import UnitsConverter
from ..collectAndCalcHelpers import (DimensionsHelper,\
                                     addLastThetaSlice,\
                                     collectiveCollect,\
                                     get2DMesh,\
                                     slicesToIndices,\
                                     polAvg)
import numpy as np

#{{{collectAndCalcFields2D
def collectAndCalcFields2D(paths                     ,\
                           varName                   ,\
                           xSlice                    ,\
                           ySlice                    ,\
                           zSlice                    ,\
                           tSlice            = None  ,\
                           convertToPhysical = True  ,\
                           fluct             = False ,\
                           uc                = None  ,\
                           dh                = None  ,\
                           mode              = "perp",\
                           xguards           = False ,\
                           yguards           = False ,\
                           ):
    #{{{docstring
    """
    Function which calculates the time traces, and gives a 4d output

    Parameters
    ----------
    paths : tuple of strings
        The paths to collect from
    varName : str
        Variable to collect
    xSlice : [Slice|int]
        The slice of the rho if the data is to be sliced.
        An integer if "mode" is set to "pol".
    ySlice : [Slice|int]
        The slice of the z if the data is to be sliced.
        An integer if "mode" is set to "perp".
    zSlice : [Slice|int]
        The slice of the z if the data is to be sliced.
        An integer if "mode" is set to "par".
    convertToPhysical : bool
        Whether or not to convert to physical units.
    fluct : bool
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    tSlice : [None|Slice]
        Whether or not to slice the time trace
    uc : [None|UnitsConverter]
        If not given, the function will create the instance itself.
        However, there is a possibility to supply this in order to
        reduce overhead.
    dh : [None|DimensionsHelper]
        If not given, the function will create the instance itself.
        However, there is a possibility to supply this in order to
        reduce overhead.
    mode : ["perp"|"par"|"pol"|]
        * "perp" - The output field is sliced along a specific z value
        * "par"  - The output field is sliced along a specific theta value
        * "pol"  - The output field is sliced along a specific rho value
    xguards : bool
        If the ghost points in x should be collected
    xguards : bool
        If the ghost points in y should be collected

    Returns
    -------
    field2D : dict
        Dictionary with the keys:
            * varName    - A 2D array of the collected variable
            * varNamePPi - The field at pi away from the varName field
                           (Only if mode is "par")
            * "X"        - The cartesian x mesh to the field
            * "Y"        - The cartesian Y mesh to the field
            * "time"     - The time trace
            * pos        - The position of the fixed index
    """
    #}}}

    # Initialize output
    field2D = {}

    if uc is None:
        # Create the units convertor object
        uc = UnitsConverter(paths[0], convertToPhysical)
    # Toggle convertToPhysical in case of errors
    convertToPhysical = uc.convertToPhysical

    if dh is None:
        # Create the dimensions helper object
        dh = DimensionsHelper(paths[0], uc)

    # Obtain the mesh
    if mode == "perp":
        # X_RT, Y_RT
        X, Y = get2DMesh(rho = dh.rho, theta = dh.theta, mode = "RT")
    elif mode == "par":
        # X_RZ, Y_RZ
        X, Y = get2DMesh(rho = dh.rho, z = dh.z, mode = "RZ")
    elif mode == "pol":
        # X_ZT, Y_ZT
        X, Y = get2DMesh(theta = dh.theta, z = dh.z, mode = "ZT")
    else:
        message = "mode '{}' not implemented".format(mode)
        raise NotImplementedError(message)

    # Convert to indices
    xInd = slicesToIndices(paths[0], xSlice, "x", xguards=xguards)
    yInd = slicesToIndices(paths[0], ySlice, "y", yguards=yguards)
    zInd = slicesToIndices(paths[0], zSlice, "z")
    tInd = slicesToIndices(paths[0], tSlice, "t")

    collectGhost = True if (xguards or yguards) else False

    # Collect
    if not(fluct):
        varTime = collectiveCollect(\
                        paths                      ,\
                        (varName, "t_array")       ,\
                        collectGhost = collectGhost,\
                        tInd         = tInd        ,\
                        yInd         = yInd        ,\
                        xInd         = xInd        ,\
                        zInd         = zInd        ,\
                        )
        var  = varTime[varName]
        time = varTime["t_array"]

        # Collect the negative
        if mode == "par":
            # Then theta index corresponding to pi
            piInd = round(var.shape[3]/2)

            if zInd > piInd:
                zNeg = zInd - piInd
            else:
                zNeg = zInd + piInd

            varPPi = collectiveCollect(\
                            paths                      ,\
                            (varName, )                ,\
                            collectGhost = collectGhost,\
                            tInd         = tInd        ,\
                            yInd         = yInd        ,\
                            xInd         = xInd        ,\
                            zInd         = zNeg        ,\
                            )
            varPPi = varPPi[varName]
    else:
        varTime = collectiveCollect(\
                            paths                      ,\
                            (varName, "t_array")       ,\
                            collectGhost = collectGhost,\
                            tInd         = tInd        ,\
                            yInd         = yInd        ,\
                            xInd         = xInd        ,\
                            zInd         = None        ,\
                            )
        var  = varTime[varName]
        time = varTime["t_array"]

    if xguards:
        # Remove the inner ghost points from the variable
        var = np.delete(var, (0), axis=1)

    # Slice in t
    if tSlice is not None:
        if tSlice.step is not None:
            var = var[::tSlice.step]

    if fluct:
        avg = polAvg(var)
        var = var - avg
        if mode == "par":
            # The negative must have the same average, but not the same
            # fluctuations
            varPPi = var - avg

    if mode == "perp":
        # Add the last theta slice
        var = addLastThetaSlice(var, var.shape[0])

    if convertToPhysical:
        var  = uc.physicalConversion(var , varName)
        time = uc.physicalConversion(time, "t")

    # Store the fields
    field2D["X"   ] = X
    field2D["Y"   ] = Y
    field2D["time"] = time

    if "pol" in mode:
        field2D[varName ] = var[:, 0, :, :]
        field2D["rhoPos"] = dh.rho[xInd[0]]
    if "perp" in mode:
        field2D[varName] = var[:, :, 0, :]
        field2D["zPos" ] = dh.z[yInd[0]]
    if "par" in mode:
        field2D[varName      ] = var   [:, :, :, 0]
        field2D[varName+"PPi"] = varPPi[:, :, :, 0]
        field2D["thetaPos"   ] = dh.theta[zInd[0]]

    return field2D
#}}}
