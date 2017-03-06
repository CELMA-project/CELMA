#!/usr/bin/env python

"""
Contains functions which help with the 3D mesh
"""

import numpy as np

#{{{addLastThetaSlice
def addLastThetaSlice(field, nFrames):
    """
    Adds the values in theta = 0 in the new point theta=2*pi
    """

    if len(field.shape) == 3:
        # Field is of x,y,z
        # Append one new dimension in the front
        field4D = np.empty((nFrames,\
                            field.shape[0],\
                            field.shape[1],\
                            field.shape[2]))
        field4D[:] = field
        field = field4D

    # Determine sizes
    oldSize     = np.array(field.shape)
    newSize     = oldSize.copy()
    newSize[-1] = oldSize[-1]+1

    # Create the new field
    newField    = np.empty(newSize)

    # Fill the new field with the old data
    # (NOTE: End-point index does not start counting on 0)
    newField[:, :, :, :oldSize[-1]] = field
    # And add the values from point theta = 0 to theta = 2*pi
    # (-1 since start count on 0)
    newField[:, :, :, newSize[-1]-1] = field[:,:,:,0]

    return newField
#}}}

#{{{get2DMesh
def get2DMesh(rho=None, thetaRad=None, z=None, mode="RT", xguards=False):
    #{{{docstring
    """
    Provides rho-theta, rho-z and theta-z meshes from provided input.

    Parameters
    ----------
    rho : [None|array]
        The rho coordinate.
    theta : [None|array]
        The theta coordinate.
        NOTE: The 2*pi point will be added before making the meshes
    z : [None|array]
        The z coordinate.
    xguards : bool
        Whether or not the xguards are present.
    mode : ["RT"|"RZ"|"ZT"]
        What mode to run:
            * RT - The rho-theta mesh will be the output
            * RZ - The rho-z mesh will be the output
            * ZT - The z-theta mesh will be the output

    Returns
    -------
    X_RT : array-2D
        Given if mode is "RT".
        The Cartesian x-mesh from the rho-theta transformation
    Y_RT : array-2D
        Given if mode is "RT".
        The Cartesian y-mesh from the rho-theta transformation
    X_RZ : array-2D
        Given if mode is "RZ".
        The Cartesian x-mesh from the rho-Z transformation (the positive
        theta position)
    Y_RZ : array-2D
        Given if mode is "RZ".
        The Cartesian y-mesh from the rho-z transformation
    X_ZT : array-2D
        Given if mode is "ZT".
        The Cartesian x-mesh from the z-theta transformation
    Y_ZT : array-2D
        Given if mode is "ZT".
        The Cartesian y-mesh from the z-theta transformation
    """
    #}}}

    if "R" in mode and xguards:
        # Remove the first point if xguards is set
        rho = rho[1:]

    if "T" in mode:
        # Add the last slice to theta
        thetaRad = np.append(thetaRad, 2.0*np.pi)

    if mode == "RT":
        # For the rho, theta plane
        THETA_RT, RHO_RT = np.meshgrid(thetaRad, rho)

        # Convert RHO and THETA to X_RT and Y_RT
        X_RT = RHO_RT*np.cos(THETA_RT)
        Y_RT = RHO_RT*np.sin(THETA_RT)

        return X_RT, Y_RT

    if mode == "RZ":
        # For the rho, z plane
        Z_RZ, RHO_RZ = np.meshgrid(z, rho)

        # Convert RHO and Z to X_RZ and Y_RZ
        X_RZ     =   RHO_RZ
        Y_RZ     =   Z_RZ

        return X_RZ, Y_RZ

    if mode == "ZT":
        # For the theta, z plane
        Z_ZT, THETA_ZT = np.meshgrid(z, thetaRad)

        # Set Z and THETA to X_ZT and Y_ZT
        X_ZT = THETA_ZT
        Y_ZT = Z_ZT

        return X_ZT, Y_ZT
#}}}
