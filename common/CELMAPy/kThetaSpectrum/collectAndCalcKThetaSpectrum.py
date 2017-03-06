#!/usr/bin/env python

"""
Contains class to calculate the magnitude spectrum of a time trace of a
spatial FFT.
"""

from ..growthRates import CollectAndCalcGrowthRates
from ..superClasses import CollectAndCalcPointsSuperClass
from ..collectAndCalcHelpers import timeAvg
import numpy as np




# FIXME: THIS IS BEING REWRITTEN, DO NOT EXPECT TO MAKE ANY SENSE OF THIS




#{{{CollectAndCalcKThetaSpectrum
class CollectAndCalcKThetaSpectrum(CollectAndCalcPointsSuperClass):
    """
    Class for collecting and calculating the magnitude spectrum
    """

    #{{{constructor
    def __init__(self             ,\
                 collectPaths     ,\
                 varName          ,\
                 convertToPhysical,\
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
        indicesArgs : tuple
            Tuple containing the indices.
        indicesKwargs : dict
            Keyword arguments to use when setting the indices for
            collecting.
        """
        #}}}

        # Set the member data
        self._collectPaths      = collectPaths
        self._varName           = varName
        self._convertToPhysical = convertToPhysical
        self._indicesArgs       = indicesArgs
        self._indicesKwargs     = indicesKwargs
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Returns
        -------
        kThetaSpectrum : dict
            The dictionary containing the mode spectra.
            The keys of the dict is on the form (rho,z), with a new dict
            as the value. This dict contains the following keys:
                * "modeAvg" - The average magnitude of the mode
                * "modeStd" - The standard deviation of the magnitude of
                              the mode
                * "modeNr"  - The mode number
        uc : Units Converter
            The units converter used when obtaining the fourier modes.
            Needed in the plotting routine.
        """
        #}}}

        steadyStatePath = self._indicesKwargs["steadyStatePath"]
        fm, _, uc =\
                CollectAndCalcGrowthRates.collectAndCalcFourierModes(\
                                   self._collectPaths         ,\
                                   self._varName              ,\
                                   self._convertToPhysical    ,\
                                   steadyStatePath            ,\
                                   self._indicesArgs          ,\
                                   self._indicesKwargs        ,\
                                   calcAngularFrequency = False,\
                                      )

        # Initialize the output
        kThetaSpectrum = {}

        for key in fm.keys():
            # Pop variables not needed
            fm[key].pop(self._varName)
            fm[key].pop("time")

            # Remove the offset mode
            fm[key][self._varName + "Magnitude"] =\
                    fm[key][self._varName + "Magnitude"][:, 1:]

            modes    = fm[key][self._varName + "Magnitude"]
            nModes   = modes.shape[1]
            modesAvg = np.empty(nModes)
            modesStd = np.empty(nModes)
            mNr      = np.empty(nModes)
            for i in range(nModes):
                # Bloat the modes so that timeAvg can read them
                mode    = np.expand_dims(\
                            np.expand_dims(\
                                np.expand_dims(modes[:,i],axis=-1),\
                                axis=-1),\
                            axis=-1)
                modeAvg = timeAvg(mode)
                modeStd = np.sqrt(timeAvg((mode - modeAvg)**2.0))

                # Flatten again after timeAvg
                modesAvg[i] = modeAvg.flatten()
                modesStd[i] = modeStd.flatten()
                mNr     [i] = i+1

            fm[key].pop(self._varName + "Magnitude")

            # Rename to mSpec
            kThetaSpectrum[key] = {\
                                      "modeAvg" : modesAvg,\
                                      "modeStd" : modesStd,\
                                      "modeNr"  : mNr     ,\
                                     }

        return kThetaSpectrum, uc
    #}}}
#}}}
