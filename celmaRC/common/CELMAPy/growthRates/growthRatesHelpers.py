#!/usr/bin/env python

"""
Contains helper function for the growth rates
"""

from ..plotHelpers import PlotHelper
import numpy as np
import matplotlib.pylab as plt
from matplotlib.gridspec import GridSpec
import os

#{{{linRegOfExp
def linRegOfExp(x,y):
    #{{{docstring
    """
    Calculates the gradient of an exponential function using linear
    regression.

    NOTE: There seem to be a lot of confusion of how to estimate
          uncertainties. To be consise, we will here follow
          "An introduction to error analysis" by Tylor, J.R, and
          not use what is used in polyfit or stats.linregress.

    NOTE: We are fitting an exponential function. The uncertainties in the
          gradient is strictly speaking NOT given in (8.17).
          We should use the weigthed fit, using (8.34) and the
          equation of p 199 to  estimate the uncertainties in
          sigma_b.
          However, we are here only interested to get the spread in the
          data of the straigth line (the standard deviation of the
          straigth line), and equation (8.17) will be used. If we were
          to calculate the spread in the exponential data (assuming that
          the uncertainties were constant), we would have

            >>> # Calculation of weigths
            >>> # We have that
            >>> # lnY = ln(y).
            >>> # We have assumed that the uncertainties in y is equally
            >>> # uncertain. From error propagation, we then get.
            >>> sigmaY = sigmaLnY/y
            >>> w      = 1/sigmaY**2
            >>> # Regression of exponential
            >>> DeltaExp = sum(w)*sum(w*x**2) - (sum(w*x))**2
            >>> BExp     = (sum(w)*sum(w*lnY*x) - sum(w*lnY)*sum(w*x))/DeltaExp
            >>> sigmaB   = np.sqrt(sum(w)/DeltaExp)

          Note also that there is not a huge difference between B and BExp

    Parameters
    ----------
    x : array
        The x data. Should be without measuring uncertainties.
    y : array
        The y data. It will be assumed that the uncertainties in y is
        equally uncertain. This means that the uncertainties in ln(y)
        will scale as sigmaLnY/y

    Returns
    -------
    B : float
        The exponential growth rate
    sigmaB : float
        The uncertainties in BExp (see notes above)
    """
    #}}}

    # Take the logarithm
    lnY = np.log(y)
    # Calculation of sigmaLnY
    N        = len(x)
    Delta    = N*sum(x**2) - (sum(x))**2
    A        = (sum(x**2)*sum(lnY) - sum(x)*sum(x*lnY))/Delta
    B        = (N*sum(lnY*x) - sum(x)*sum(lnY))/Delta
    sigmaLnY = np.sqrt((1/(N-2))*sum((lnY - A - B*x)**2))
    sigmaB   = sigmaLnY*np.sqrt(N/Delta)

    return B, sigmaB
#}}}
