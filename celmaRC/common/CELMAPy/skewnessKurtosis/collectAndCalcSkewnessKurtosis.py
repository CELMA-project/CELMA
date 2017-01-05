#!/usr/bin/env python

"""
Contains class for collecting and calculating the skewness and kurtosis
of the time traces.
"""

from ..timeTrace import CollectAndCalcTimeTrace
from scipy.stats import kurtosis, skew

#{{{CollectAndCalcSkewnessKurtosis
class CollectAndCalcSkewnessKurtosis(CollectAndCalcTimeTrace):
    """
    Class for collecting and calcuating skewness and kurtosis of the time traces

    Skewness
    --------
    Negative skew: The left tail is longer; the mass of the distribution
    is concentrated on the right of the figure.
    For a pure Gaussian the skewness is 0.

    Kurtosis
    --------
    The kurtosis of any univariate normal distribution is 3.
    Kurtosis less than 3: Platykurtic: The distribution produces fewer
    and less extreme outliers than a gaussian.
    Kurtosis greater than 3: Leptokurtic: Tails approaches zero slower than a
    Gaussian (more extreme outliners).
    """

    #{{{constructor
    def __init__(self   ,\
                 *args  ,\
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
    #{{{calcSkewnessKurtosis
    def calcSkewnessKurtosis(timeTraces):
        #{{{docstring
        """
        Function which calculates the skewness and kurtosis.

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
        skewKurt : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varNameSkew:skewness, varNameKurt:kurtosis}
            NOTE: Fisher's kurtosis (excess) is used, where 3 is subtracted from
                  Pearson's definition
        """
        #}}}

        # Initialize the output
        skewKurt = {}

        # Obtain the varName
        ind  = tuple(timeTraces.keys())[0]
        keys = timeTraces[ind].keys()
        varName = tuple(var for var in keys if var != "time")[0]

        # Make the keys
        skewKey = "{}Skew".format(varName)
        kurtKey = "{}Kurt".format(varName)

        # Obtain the skewness and kurtosis
        for key in timeTraces.keys():
            # Initialize the dict
            skewKurt[key] = {}

            # Calculate the skewness and kurtosis
            # NOTE: Default for kurtosis is to have Fisher = True and
            #       bias correction on
            skewKurt[key][skewKey] = skew(timeTraces[key][varName])
            skewKurt[key][kurtKey] = kurtosis(timeTraces[key][varName])

        return skewKurt
    #}}}
#}}}
