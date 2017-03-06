#!/usr/bin/env python

"""
Contains super class for setting data for collection.
"""

from ..collectAndCalcHelpers import DimensionsHelper
from ..unitsConverter import UnitsConverter

#{{{CollectAndCalcSuperClass
class CollectAndCalcSuperClass(object):
    """
    The parent collect and calc class
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
        Super constructor for collection and calculatio.

        This constructor will:
            * Set the member data
            * Create the UnitsConverter
            * Create the DimensionsHelper.

        Parameters
        ----------
        collectPaths : tuple of strings
            The paths to collect from
        convertToPhysical : bool
            Whether or not to convert to physical units.
        xguards : bool
            If the ghost points in x should be collected.
        xguards : bool
            If the ghost points in y should be collected.
        uc : [None|UnitsConverter]
            If not given, the function will create the instance itself.
            However, there is a possibility to supply this in order to
            reduce overhead.
        dh : [None|DimensionsHelper]
            If not given, the function will create the instance itself.
            However, there is a possibility to supply this in order to
            reduce overhead.
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

        self._notCalled = ["setVarName"]
    #}}}

    #{{{setVarName
    def setVarName(self, varName):
        #{{{docstring
        """
        Sets the varName.

        Parameters
        ----------
        varName : str
            Variable to collect.
        """
        #}}}
        if "setVarName" in self._notCalled:
            self._notCalled.remove("setVarName")
        self._varName = varName
    #}}}
#}}}
