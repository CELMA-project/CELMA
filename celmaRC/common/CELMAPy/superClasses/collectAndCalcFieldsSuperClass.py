#!/usr/bin/env python

"""
Contains super class for setting data for 1D and 2D fields collection.
"""

from .collectAndCalcSuperClass import CollectAndCalcSuperClass

#{{{CollectAndCalcFieldsSuperClass
class CollectAndCalcFieldsSuperClass(CollectAndCalcSuperClass):
    """
    Super class for field collection.

    Handles setting of data for field collection.
    """

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        Super constructor for collection and calculation of fields.

        This constructor will:
            * Call the parent class

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._notCalled.append("setSlices")
    #}}}

    #{{{setSlice
    def setSlice(self, xSlice, ySlice, zSlice, tSlice):
        #{{{docstring
        """
        Sets the slices

        Parameters
        ----------
        xSlice : [None|slice|int]
            The slice of the rho if the data is to be sliced.
            If int, a constant slice will be used.
        ySlice : [None|slice|int]
            The slice of the z if the data is to be sliced.
            If int, a constant slice will be used.
        zSlice : [None|slice|int]
            The slice of the z if the data is to be sliced.
            If int, a constant slice will be used.
        tSlice : [None|Slice]
            Whether or not to slice the time trace
        """
        #}}}
        if "setSlices" in self._notCalled:
            self._notCalled.remove("setSlices")

        self._xSlice = xSlice
        self._ySlice = ySlice
        self._zSlice = zSlice
        self._tSlice = tSlice
    #}}}
#}}}
