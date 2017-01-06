#!/usr/bin/env python

"""
Contains class for collecting and calculating the radial profiles
"""

from ..fields1D import CollectAndCalcFields1D
from ..collectAndCalcHelpers import (calcN, calcUIPar, calcUEPar,\
                                     polAvg, timeAvg)
import numpy as np

#{{{CollectAndCalcRadialProfile
class CollectAndCalcRadialProfile(object):
    """
    Class for collecting and calcuating radial profiles
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
        zInd = 0
        self._slices = (None, yInd, zInd, tSlice)

        # Placeholder for uc and dh
        self.uc = None
        self.dh = None
    #}}}

    #{{{collectWrapper
    def collectWrapper(self, paths, varName):
        #{{{docstring
        """
        Wrapper function around the collect routine.

        Sets the units converter if unset

        Parameters
        ----------
        paths : tuple
            Tuple of the collect paths
        varName : str
            Name of the variable to collect

        Returns
        -------
        dict1D : dict
            Dictionary with the added variable.
            For details, see the documentation of the CollectAndCalcFields1D
            class.
        """
        #}}}
        ccf1D = CollectAndCalcFields1D(\
                paths,\
                mode = "radial",\
                processing = None,\
                return2d   = False,\
                convertToPhysical = self._convertToPhysical)

        # Set the units converter if not set
        if self.uc is None:
            self.uc = ccf1D.uc
        # Set the dimension helper if not set
        if self.dh is None:
            self.dh = ccf1D.dh

        ccf1D.setSlice(*self._slices)
        specialCollects = ("n", "uIPar", "uEPar")
        if not(varName in specialCollects):
            ccf1D.setVarName(varName)
            dict1D = ccf1D.executeCollectAndCalc()
        else:
            dict1D = self._specialCollect(varName, ccf1D)

        return dict1D
    #}}}

    #{{{_specialCollect
    def _specialCollect(self, varName, ccf1D):
        #{{{docstring
        """
        Calculate non-collectable variables

        Parameters
        ----------
        varName : str
            Name of the variable to calculate
        ccf1D : CollectAndCalcFields1D
            Object to perform helper collections

        Returns
        -------
        dict1D : dict
            Dictionary with the added variable.
            For details, see the documentation of the CollectAndCalcFields1D
            class.
        """
        #}}}

        ccf1D.setVarName("lnN")
        dict1D = ccf1D.executeCollectAndCalc()
        dict1D.update({"n" : calcN(dict1D["lnN"],\
                                   not(ccf1D.convertToPhysical),\
                                   ccf1D.uc)})
        if varName != "n":
            ccf1D.setVarName("momDensPar")
            dict1D.update(ccf1D.executeCollectAndCalc())
            dict1D.update({"uIPar":\
                           calcUIPar(dict1D["momDensPar"],dict1D["n"])})
        if varName == "uEPar":
            ccf1D.setVarName("jPar")
            dict1D.update(ccf1D.executeCollectAndCalc())
            dict1D.update({"uEPar": calcUEPar(dict1D["uIPar"],\
                                              dict1D["jPar"],\
                                              dict1D["n"],\
                                              not(ccf1D.convertToPhysical))})
        return dict1D
    #}}}

    @staticmethod
    #{{{calcAvgFluctStd
    def calcAvgFluctStd(var):
        #{{{docstring
        """
        Calculates the averge, fluctuation and standard deviation

        Parameters
        ----------
        var : array-4d
            Array to use in the calcualtions

        Returns
        -------
        avg : array-4d
            The average of var
        fluct : array-4d
            The fluct of var
        std : array-4d
            The standar deviation
        """
        #}}}
        # Cacluclate the average
        avg = polAvg(timeAvg(var))
        # Caclculate the fluctuation
        fluct = var - avg
        # Calculate the standard deviation
        # This follows from the definition of the standard deviation
        std = np.sqrt(polAvg(timeAvg(fluct**2.0)))

        return avg, fluct, std
    #}}}
#}}}
