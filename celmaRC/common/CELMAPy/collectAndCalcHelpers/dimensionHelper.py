#!/usr/bin/env python

""" Contains the DimensionsHelper class """

from .gridSizes import getUniformSpacing, getMYG, getMXG, getSizes
import numpy as np

#{{{DimensionsHelper
class DimensionsHelper(object):
    """Contains routines for dealing with the dimensions."""

    #{{{__init__
    def  __init__(self                     ,\
                  path                     ,\
                  unitsConverter           ,\
                  xguards           = False,\
                  yguards           = False):
        #{{{docstring
        """
        The constructor for DimensionsHelper, which:

        * Collect the coordinates
        * Convert the dimensions and the spacings

        Parameters
        ----------
        path : str
            The path to collect from.
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        unitsConverter : UnitsConverter
            UnitsConverter object which contains the conversion factors.
        """
        #}}}

        # Set the member data
        self._path              = path
        self._unitsConverter    = unitsConverter
        self._xguards           = xguards
        self._yguards           = yguards
        self._convertToPhysical = unitsConverter.convertToPhysical

        # Get the coordinates
        self._getCoordinates()

        # Convert the coordinates
        self._convertCoordinates()
    #}}}

    #{{{_getCoordinates
    def _getCoordinates(self):
        #{{{docstring
        """
        Obtains and calculates the coordinates.

        The coordinates are available through the CoordinatesHelper
        object, and consists of
            * self.rho   - The rho coordinate
            * self.dx    - The grid spacing in rho
            * self.MXG   - The number of ghost points in rho
            * self.z     - NOTE THE NAME! The z coordinate
            * self.dy    - NOTE THE NAME! The grid spacing in z
            * self.MYG   - The number of ghost points in z
            * self.theta - The theta coordinate
            * self.dz    - NOTE THE NAME! The grid spacing in theta
        """
        #}}}

        #{{{rho
        dx = getUniformSpacing(self._path           ,\
                               "x"                  ,\
                               xguards=self._xguards,\
                               yguards=self._yguards)

        MXG = getMXG(self._path)

        nPoints = dx.shape[0]
        dx      = dx[0,0]

        if self._xguards:
            innerPoints = nPoints - 2*MXG
        else:
            innerPoints = nPoints

        rho = dx * np.array(np.arange(0.5, innerPoints))

        if self._xguards:
            # Insert the first and last grid point
            rho = np.insert(rho, 0, - 0.5*dx)
            rho = np.append(rho, rho[-1] + dx)

        self.rho = rho
        self.dx  = dx
        self.MXG = MXG
        #}}}

        #{{{z
        dy = getUniformSpacing(self._path           ,\
                               "y"                  ,\
                               xguards=self._xguards,\
                               yguards=self._yguards)

        MYG = getMYG(self._path)

        nPoints  = dy.shape[1]
        dy = dy[0,0]

        if self._yguards:
            innerPoints = nPoints - 2*MYG
        else:
            innerPoints = nPoints

        z = dy * np.array(np.arange(0.5, innerPoints))

        if self._yguards:
            # Insert the first and last grid point
            z = np.insert(z, 0, - 0.5*dy)
            z = np.append(z, z[-1] + dy)

        self.z   = z
        self.dy  = dy
        self.MYG = MYG
        #}}}

        #{{{theta
        dz = getUniformSpacing(self._path, "z")

        innerPoints = getSizes(self._path, "z")

        theta = dz * np.array(np.arange(0.0, innerPoints))

        # Convert to degrees
        self.thetaRad = theta
        theta = theta * (180/np.pi)

        self.thetaGrad = theta
        self.dz        = dz
        #}}}
    #}}}

    #{{{_convertCoordinates
    def _convertCoordinates(self):
        #{{{docstring
        """
        Converts the following units if convertToPhysical is True:
            * self.t
            * self.dt
            * self.rho
            * self.dx
            * self.z
            * self.dy
        """
        #}}}

        if self._convertToPhysical:
            self.rho = self._unitsConverter.physicalConversion(self.rho, "rho")
            self.dx  = self._unitsConverter.physicalConversion(self.dx , "rho")
            self.z   = self._unitsConverter.physicalConversion(self.z  , "z"  )
            self.dy  = self._unitsConverter.physicalConversion(self.dy , "z"  )
    #}}}
#}}}
