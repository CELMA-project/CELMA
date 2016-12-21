#!/usr/bin/env python

"""
Contains the super class driver for fields 1D and fields 2D
"""

from .driverPostProcessingSuperClass import DriverPostProcessingSuperClass

#{{{DriverPlotFieldsSuperClass
class DriverPlotFieldsSuperClass(DriverPostProcessingSuperClass):
    """
    The parent driver of 1D and 2D field plotting
    """

    #{{{Constructor
    def __init__(self           ,\
                 *args          ,\
                 xguards = False,\
                 yguards = False,\
                 xSlice  = None ,\
                 ySlice  = None ,\
                 zSlice  = None ,\
                 tSlice  = None ,\
                 xInd    = None ,\
                 yInd    = None ,\
                 zInd    = None ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Sets the common memberdata

        Parameters
        ----------
        *args : str
            See parent class for details.
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
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

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
    #}}}
#}}}
