#!/usr/bin/env python

"""
Init-file for plotting
"""

from .driversCombined import combinedDriver, combined1D2D
from .drivers1D import single1DDriver, parDriver, perpDriver, parPerpDriver
from .drivers2D import (single2DDriver       ,\
                        allMainFields2DDriver,\
                        lnN2DDriver          ,\
                        uEPar2DDriver        ,\
                        uIPar2DDriver        ,\
                        vortD2DDriver        ,\
                        vort2DDriver         ,\
                        phi2DDriver          ,\
                       )

import matplotlib.pyplot as plt

# Set proper backend
try:
    plt.figure(0)
except RuntimeError:
    plt.switch_backend('Agg')
    plt.figure(0)
plt.close(0)
