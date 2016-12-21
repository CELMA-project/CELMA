#!/usr/bin/env python

"""
Contains class for collecting and calculating the 2D fields
"""

from ..collectAndCalcHelpers import DimensionsHelper
from ..unitsConverter import UnitsConverter

#{{{CollectAndCalcFieldsSuperClass
class CollectAndCalcFieldsSuperClass(object):
    """
    Super class for field collection.

    Handles common plot options and saving.
    """

    #{{{constructor
    def __init__(self                      ,\
                 collectPaths              ,\
                 convertToPhysical = True  ,\
                 xguards           = False ,\
                 yguards           = False ,\
                 uc                = None  ,\
                 dh                = None  ,\
                ):
        #{{{docstring
        """
        Super constructor for collection and calculation of fields.

        This constructor will:
            * Set the member data
            * Create the UnitsConverter
            * Create the DimensionsHelper.

        Parameters
        ----------
        collectPaths : tuple of strings
            The paths to collect from
        varName : str
            Variable to collect
        convertToPhysical : bool
            Whether or not to convert to physical units.
        uc : [None|UnitsConverter]
            If not given, the function will create the instance itself.
            However, there is a possibility to supply this in order to
            reduce overhead.
        dh : [None|DimensionsHelper]
            If not given, the function will create the instance itself.
            However, there is a possibility to supply this in order to
            reduce overhead.
        xguards : bool
            If the ghost points in x should be collected
        xguards : bool
            If the ghost points in y should be collected
        """
        #}}}

        self._collectPaths = collectPaths

        self._xguards = xguards
        self._yguards = yguards

        if uc is None:
            # Create the units convertor object
            uc = UnitsConverter(collectPaths[0], convertToPhysical)
        # Toggle convertToPhysical in case of errors
        self.convertToPhysical = uc.convertToPhysical

        if dh is None:
            # Create the dimensions helper object
            dh = DimensionsHelper(collectPaths[0], uc)

        self.uc = uc
        self._dh = dh

        self._notCalled = ["setVarName", "setSlices"]
    #}}}

    #{{{setVarName
    def setVarName(self, varName):
        #{{{docstring
        """
        Sets the varName

        Parameters
        ----------
        varName : str
            Variable to collect
        """
        #}}}
        if "setVarName" in self._notCalled:
            self._notCalled.remove("setVarName")
        self._varName = varName
    #}}}

    #{{{setSlice
    def setSlice(self, xSlice, ySlice, zSlice, tSlice):
        #{{{docstring
        """
        Sets the slices

        Parameters
        ----------
        xSlice : [Slice|int]
            The slice of the rho if the data is to be sliced.
            If int, a constant slice will be used.
        ySlice : [Slice|int]
            The slice of the z if the data is to be sliced.
            If int, a constant slice will be used.
        zSlice : [Slice|int]
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
