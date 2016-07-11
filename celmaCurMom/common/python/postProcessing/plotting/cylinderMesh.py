#!/usr/bin/env python

"""
Class which creates a the points in a cylinder.
"""

import numpy as np

#{{{class CylinderMesh
class CylinderMesh(object):
    """
    Provides the grid of a cylinder

    Also provides a method to fill the last theta slice.
    """

    #{{{Constructor
    def __init__(self   ,\
                 rho    ,\
                 theta  ,\
                 z      ,\
                 xguards,\
                ):
        """
        Provides rho-theta and rho-z meshes from provided input.

        theta will be added one additional point to be consistent with
        addLastThetaSlice
        """

        # Remove the first point if xguards is set
        if xguards:
            rho = rho[1:]

        # Add the last slice to theta
        theta = np.append(theta, 2.0*np.pi)

        # For the rho, theta plane
        THETA_RT, RHO_RT = np.meshgrid(theta, rho)
        # For the rho, z plane
        Z_RZ, RHO_RZ = np.meshgrid(z, rho)

        # Convert RHO and THETA to X_RT and Y_RT
        self.X_RT = RHO_RT*np.cos(THETA_RT)
        self.Y_RT = RHO_RT*np.sin(THETA_RT)

        # Convert RHO and Z to X_RZ and Y_RZ
        self.X_RZ     =   RHO_RZ
        self.Y_RZ     =   Z_RZ
        self.X_RZ_NEG = - RHO_RZ
    #}}}

    #{{{addLastThetaSlice
    def addLastThetaSlice(self, field, nFrames):
        """
        Adds the values in theta = 0 in the new point theta=2*pi
        """

        if len(field.shape) == 3:
            # Field is of x,y,z
            # Append one new dimension in the front
            field4D = np.empty((nFrames ,field.shape[0], field.shape[1], field.shape[2]))
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
#}}}
