#!/usr/bin/env python

"""
Contains class for collecting and calculating the probability density
function of time traces.
"""

from ..timeTrace import CollectAndCalcTimeTrace
import numpy as np

#{{{CollectAndCalcPDF
class CollectAndCalcPDF(CollectAndCalcTimeTrace):
    """
    Class for collecting and calcuating PDF of the time traces

    Probability density function
    ----------------------------
    Probability that the measurement falls within an infinite small
    interval.
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
    #{{{calcPDF
    def calcPDF(timeTraces):
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
        PDF : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varPDFX:pdfX, varPDFY:pdfY}
        """
        #}}}

        # Initialize the output
        PDF = {}

        # Obtain the varName
        ind  = tuple(timeTraces.keys())[0]
        keys = timeTraces[ind].keys()
        varName = tuple(var for var in keys if var != "time")[0]

        # Make the keys
        xKey = "{}PDFX".format(varName)
        yKey = "{}PDFY".format(varName)

        # Histogram counts the occurences of values within a specific interval
        # Density normalizes so that the integral (of the continuous variable)
        # equals one, note that the sum of histograms is not necessarily 1)
        # http://docs.scipy.org/doc/numpy/reference/generated/numpy.histogram.html
        # http://stackoverflow.com/questions/36150257/probability-distribution-function-python/36248810
        for key in timeTraces.keys():
            # Initialize the PDF
            PDF[key] = {}

            # Calculate pdfY
            PDF[key][yKey], bins =\
                np.histogram(timeTraces[key][varName],\
                             bins="auto",\
                             density=True)

            # Initialize x
            PDF[key][xKey] = np.zeros(PDF[key][yKey].size)

            # Only the bin edges are saved. Interpolate to bin center
            for k in range(PDF[key][yKey].size):
                PDF[key][xKey][k] = 0.5*(bins[k]+bins[k+1])

        # NOTE: If timeTraces was converted to physical units, then PDF is
        #       in physical units as well
        return PDF
    #}}}
#}}}
