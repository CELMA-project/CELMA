#!/usr/bin/env python

"""
Contains functions for calculating the ExB velocities
"""

from ..fields1D import CollectAndCalcFields1D
from .averages import polAvg
from .derivatives import DDZ
import numpy as np

#{{{calcRadialExBPoloidal
def calcRadialExBPoloidal(collectPaths, slices, dh,\
                          mode="fluct", convertToPhysical = True):
    #{{{docstring
    """
    Calculates the radialExB velocity in a poloidal profile

    Parameters
    ----------
    collectPaths : tuple
        Tuple from where to collect
    slices : tuple
        Tuple the indices to use.
        On the form (xInd, yInd, tSlice)
    dh : DimensionsHelper
        DimensionHelper object (used to find the Jacobian)
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

    # Calculate the radial flux
    J = np.array((dh.rho,))
    DDZPhi = DDZ(phi, J)
    if mode == "fluct":
        DDZPhi = (DDZPhi - polAvg(DDZPhi))

    radialExB = DDZPhi
    
    return radialExB, phiDict.pop("time")
#}}}
