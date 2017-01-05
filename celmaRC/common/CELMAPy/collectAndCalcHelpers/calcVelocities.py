#!/usr/bin/env python

"""
Contains functions for calculating the ExB velocities
"""

from ..fields1D import CollectAndCalcFields1D
from .averages import polAvg
from .derivatives import DDX, DDZ
import numpy as np

#{{{calcRadialExBPoloidal
def calcRadialExBPoloidal(collectPaths, slices,\
                          mode="fluct", convertToPhysical = True):
    #{{{docstring
    """
    Calculates the radial ExB velocity in a poloidal profile

    Parameters
    ----------
    collectPaths : tuple
        Tuple from where to collect
    slices : tuple
        Tuple the indices to use.
        On the form (xInd, yInd, tSlice)
    mode : ["normal"|"fluct"]
        Whether to look at fluctuations or normal data
    convertToPhysical : bool
        Whether or not to convert to physical

    Returns
    -------
    radialExB : array-4d
        The radial ExB in the fixed rho and z position
    time : array-1d
        The corresponding time
    """
    #}}}

    # Reform the slice
    xInd, yInd, tSlice = slices
    slices = (xInd, yInd, None, tSlice)

    # Collect phi
    ccf1D = CollectAndCalcFields1D(\
                collectPaths,\
                mode = "poloidal" ,\
                return2d = False  ,\
                convertToPhysical = convertToPhysical)

    ccf1D.setSlice(*slices)
    ccf1D.setVarName("phi")
    phiDict = ccf1D.executeCollectAndCalc()
    phi = phiDict.pop("phi")

    # Calculate the radial ExB
    dh = ccf1D.getDh()
    J = np.array((dh.rho,))
    DDZPhi = DDZ(phi, J)
    if mode == "fluct":
        DDZPhi = (DDZPhi - polAvg(DDZPhi))

    radialExB = DDZPhi

    return radialExB, phiDict.pop("time")
#}}}

#{{{calcPoloidalExBConstZ
def calcPoloidalExBConstZ(collectPaths, slices,\
                          mode="fluct", convertToPhysical = True):
    #{{{docstring
    """
    Calculates the poloidal ExB velocity for a parallel plane

    Parameters
    ----------
    collectPaths : tuple
        Tuple from where to collect
    slices : tuple
        Tuple the indices to use.
        On the form (yInd, tSlice)
    mode : ["normal"|"fluct"]
        Whether to look at fluctuations or normal data
    convertToPhysical : bool
        Whether or not to convert to physical

    Returns
    -------
    poloidalExB : array-4d
        The poloidal ExB for a fixed parallel plane
    time : array-1d
        The corresponding time
    """
    #}}}

    # Reform the slice
    zInd, tSlice = slices
    slices = (None, yInd, None, tSlice)

    # Collect phi
    ccf1D = CollectAndCalcFields1D(\
                collectPaths,\
                mode = "radial" ,\
                return2d = False  ,\
                convertToPhysical = convertToPhysical)

    ccf1D.setSlice(*slices)
    ccf1D.setVarName("phi")
    phiDict = ccf1D.executeCollectAndCalc()
    phi = phiDict.pop("phi")

    # Calculate the poloidal ExB
    dh = ccf1D.getDh()
    DDXPhi = DDX(phi, dh.dx)
    if mode == "fluct":
        DDXPhi = (DDXPhi - polAvg(DDXPhi))

    poloidalExB = DDXPhi

    return poloidalExB, phiDict.pop("time")
#}}}
