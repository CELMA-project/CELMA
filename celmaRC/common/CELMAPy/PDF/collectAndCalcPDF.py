#!/usr/bin/env python

"""
Contains the PDF calculation
"""

from ..timeTrace import calcTimeTrace
import numpy as np

#{{{calcPDF
def calcPDF(paths                      ,\
            varName                    ,\
            xInd                       ,\
            yInd                       ,\
            zInd                       ,\
            convertToPhysical = True   ,\
            mode              = "fluct",\
            tSlice            = None):
    #{{{docstring
    """
    Function which calculates the probability distribution function.

    Probability distribution function
    ---------------------------------
    Probability that the measurement falls within an infinite small
    interval.

    Parameters
    ----------
    paths : tuple of strings
        The paths to collect from
    varName : str
        Variable to collect
    xInd : tuple of ints
        A tuple of the xInds to collect use when collecting.
    yInd : tuple of ints
        The same as xInd, but for the y-index.
    zInd : tuple of ints
        The same as xInd, but for the z-index.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is used.
        If mode is "fluct" the fluctuations are used.
    tSlice : [None|Slice}
        Whether or not to slice the time trace

    Returns
    -------
    PDF : dict
        Dictionary where the keys are on the form "rho,theta,z".
        The value is a dict containing of
        {varPDFX:pdfX, varPDFY:"pdfY"}
    """
    #}}}

    # Call the time calcTimeTrace function in order to get a time trace
    timeTraces =\
        calcTimeTrace(paths                                ,\
                      varName                              ,\
                      xInd                                 ,\
                      yInd                                 ,\
                      zInd                                 ,\
                      convertToPhysical = convertToPhysical,\
                      mode              = mode             ,\
                      tSlice            = tSlice           ,\
                      )

    # Initialize the output
    PDF = {}

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
            np.histogram(timeTraces[key][varName], bins="auto", density=True)

        # Initialize x
        PDF[key][xKey] = np.zeros(PDF[key][yKey].size)

        # Only the bin edges are saved. Interpolate to the center of the bin
        for k in range(PDF[key][yKey].size):
            PDF[key][xKey][k] = 0.5*(bins[k]+bins[k+1])

    # NOTE: If timeTraces was converted to physical units, then PDF is
    #       in physical units as well
    return PDF
#}}}
