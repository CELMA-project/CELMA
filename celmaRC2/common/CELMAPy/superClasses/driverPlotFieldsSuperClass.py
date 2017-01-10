#!/usr/bin/env python

"""
Contains the super class driver for fields 1D and fields 2D
"""

from .driverSuperClass import DriverSuperClass

#{{{DriverPlotFieldsSuperClass
class DriverPlotFieldsSuperClass(DriverSuperClass):
    """
    The parent driver of 1D and 2D field plotting
    """

    #{{{Constructor
    def __init__(self                    ,\
                 *args                   ,\
                 convertToPhysical = True,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class

        Parameters
        ----------
        *args : positional arguments
            See parent class for details.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructors of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self.convertToPhysical = convertToPhysical
    #}}}

    #{{{setGuardSlicesAndIndices
    def setGuardSlicesAndIndices(self           ,\
                                 xguards = False,\
                                 yguards = False,\
                                 xSlice  = None ,\
                                 ySlice  = None ,\
                                 zSlice  = None ,\
                                 tSlice  = None ,\
                                 xInd    = None ,\
                                 yInd    = None ,\
                                 zInd    = None ):
        #{{{docstring
        """
        Sets the guard, the slices and the indices

        Parameters
        ----------
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        xSlice : slice
            How the data will be sliced in x.
        ySlice : slice
            How the data will be sliced in y.
        zSlice : slice
            How the data will be sliced in z.
        tSlice : slice
            How the data will be sliced in t.
        xInd : int
            Fixed rho index.
        yInd : int
            Fixed z index.
        zInd : int
            Fixed theta index.
        """
        #}}}

        # Set the member data
        self._xguards = xguards
        self._yguards = yguards
        self._xSlice  = xSlice
        self._ySlice  = ySlice
        self._zSlice  = zSlice
        self._tSlice  = tSlice
        self._xInd    = xInd
        self._yInd    = yInd
        self._zInd    = zInd
    #}}}}
#}}}
