#!/usr/bin/env python

"""
Contains class to calculate the magnitude spectrum of a time trace of a
spatial FFT.
"""

from ..growthRates import CollectAndCalcGrowthRates
from ..superClasses import CollectAndCalcPointsSuperClass

#{{{CollectAndCalcMagnitudeSpectrum
class CollectAndCalcMagnitudeSpectrum(CollectAndCalcPointsSuperClass):
    """
    Class for collecting and calculating the magnitude spectrum
    """

    #{{{constructor
    def __init__(self             ,\
                 collectPaths     ,\
                 varName          ,\
                 convertToPhysical,\
                 steadyStatePath  ,\
                 indicesArgs      ,\
                 indicesKwargs    ,\
                 ):
        #{{{docstring
        """
        This constructor will:
            * Set member data

        Parameters
        ----------
        collectPaths : tuple
            Paths to collect from.
            The corresponind 't_array' of the paths must be in ascending order.
        varName : str
            Name of variable to find the growth rates of.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        steadyStatePath : str
            String containing the steady state path
        indicesArgs : tuple
            Tuple containing the indices.
        indicesKwargs : dict
            Keyword arguments to use when setting the indices for
            collecting.
        """
        #}}}

        # Set the member data
        self._collectPaths      = collectPaths
        self._steadyStatePath   = steadyStatePath
        self._collectPaths      = collectPaths
        self._varName           = varName
        self._convertToPhysical = convertToPhysical
        self._steadyStatePath   = steadyStatePath
        self._indicesArgs       = indicesArgs
        self._indicesKwargs     = indicesKwargs
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
# FIXME:
        #{{{docstring
        """
        Returns
        -------
        fm : dict
            The dictionary containing the fourier modes, containing the
            following keys:
                * foo
        positionTuple : tuple
            The tuple containing (rho, z).
            Needed in the plotting routine.
        uc : Units Converter
            The units converter used when obtaining the fourier modes.
            Needed in the plotting routine.
        """
        #}}}

        fm, positionTuple, uc =\
                CollectAndCalcGrowthRates.collectAndCalcFourierModes(\
                                   (self._collectPaths,)      ,\
                                   self._varName              ,\
                                   self._convertToPhysical    ,\
                                   self._steadyStatePath      ,\
                                   self._indicesArgs          ,\
                                   self._indicesKwargs        ,\
                                   calcAngularVelocity = False,\
                                      )

        # Check number of modes, and how data is packed
        import pdb; pdb.set_trace()
        a = 1

        # for each mode average and find the standard deviation

        return fm, positionTuple, uc
    #}}}
#}}}
