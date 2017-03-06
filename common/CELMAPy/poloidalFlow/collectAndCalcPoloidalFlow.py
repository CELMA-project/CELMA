#!/usr/bin/env python

"""
Contains class for collecting and calculating the poloidal flows
"""

from ..collectAndCalcHelpers import DimensionsHelper
from ..calcVelocities import calcPoloidalExBConstZ
from ..unitsConverter import UnitsConverter

#{{{CollectAndCalcPoloidalFlow
class CollectAndCalcPoloidalFlow(object):
    """
    Class for collecting and calcuating poloidal flows
    """

    #{{{constructor
    def __init__(self  ,\
                 yInd  ,\
                 tSlice,\
                 convertToPhysical):
        #{{{docstring
        """
        This constructor will:
            * Set member data

        Parameters
        ----------
        yInd : int
            z-position in the cylinder
        tSlice : slice
            How the data will be sliced in time
        convertToPhysical : bool
            Whether or not to convert to physical units.
        """
        #}}}

        self._convertToPhysical = convertToPhysical
        # Notice the zInd is irrelevant, and will be not be used in the
        # collect
        self._yTSlices = (yInd, tSlice)

        # Placeholder for uc
        self.uc = None
    #}}}

    #{{{calcPoloidalExB
    def calcPoloidalExB(self, paths):
        #{{{docstring
        """
        Wrapper function around the calculation routine.

        Sets the units converter and the dimension helper

        Parameters
        ----------
        paths : tuple
            Tuple of the collect paths

        Returns
        -------
        dict1D : dict
            Dictionary with the items:
                * "uExBPoloidal" - 4d array containing the poloidalExB at
                                   the specified z-position.
                * "time"         - 1d array of the corresponding time
                * "zPos"         - float stating the fixed z position
        """
        #}}}

        poloidalExB, time =\
            calcPoloidalExBConstZ(paths         ,\
                                  self._yTSlices,\
                                  mode="normal" ,\
                                  convertToPhysical = self._convertToPhysical)

        dict1D = {"uExBPoloidal" : poloidalExB, "time":time}

        if self.uc is None:
            self.uc = UnitsConverter(paths[0], self._convertToPhysical)
            self._convertToPhysical = self.uc.convertToPhysical
            self.dh = DimensionsHelper(paths[0], self.uc)

        # Get the z position
        dict1D["zPos"] = self.dh.z[self._yTSlices[0]]

        return dict1D
    #}}}
#}}}
