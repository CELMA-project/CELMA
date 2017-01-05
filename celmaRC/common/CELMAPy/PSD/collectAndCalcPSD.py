#!/usr/bin/env python

"""
Contains class for collecting and calculating the power spectral density
function of the time traces.
"""

from ..timeTrace import CollectAndCalcTimeTrace
from scipy.signal import periodogram

#{{{CollectAndCalcPSD
class CollectAndCalcPSD(CollectAndCalcTimeTrace):
    """
    Class for collecting and calcuating PSD of the time traces

    Power spectral density
    ----------------------
        * Tells us what frequencies are present in the signal.
        * The average power of the signal is the integrated of the PSD
        * The bandwidth of the process (turbulence) is defined as
          the frequency width where the signal is within 3dB of its
          peak value.
        * PSD is a deterministic description of the spectral
          characteristic of the signal. Cannot use fourier transform on
          random variables as this is not necessarily definied etc.
        * PSD is the fourier transformed of the auto-correlation
          function, and cross spectral density is the fourier
          transformed of the cross correlation function
        * The periodogram estimate is the same as autocorrelation with a
          triangular window
    """

    #{{{constructor
    def __init__(self            ,\
                 *args           ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)
    #}}}

    @staticmethod
    #{{{calcPSD
    def calcPSD(timeTraces):
        #{{{docstring
        """
        Function which calculates the probability density function.

        Parameters
        ----------
        timeTraces : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:timeTrace, "time":time}.
            And additional key "zInd" will be given in addition to varName
            and "time" if mode is set to "fluct".
            The timeTrace is a 1d array.

        Returns
        -------
        PSD : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varPSDX:pdfX, varPSDY:pdfY}
        """
        #}}}

        # Initialize the output
        PSD = {}

        # Obtain the varName
        ind  = tuple(timeTraces.keys())[0]
        keys = timeTraces[ind].keys()
        varName = tuple(var for var in keys if var != "time")[0]

        # Make the keys
        xKey = "{}PSDX".format(varName)
        yKey = "{}PSDY".format(varName)

        # Obtain the PSD
        for key in timeTraces.keys():
            # Initialize the PSD
            PSD[key] = {}

            # Sampling frequency
            fs = 1/(timeTraces[key]["time"][1] - timeTraces[key]["time"][0])

            # window = None => window = "boxcar"
            # scaling = density gives the correct units
            PSD[key][xKey], PSD[key][yKey] =\
                periodogram(timeTraces[key][varName],\
                            fs=fs, window=None, scaling="density")

        # NOTE: If timeTraces was converted to physical units, then PSD is
        #       in physical units as well
        return PSD
    #}}}
#}}}
