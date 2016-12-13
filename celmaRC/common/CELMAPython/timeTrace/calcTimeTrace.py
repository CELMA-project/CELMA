#!/usr/bin/env python

"""
Contains classes which probes the data
"""

from .averages import polAvg
from ..plotHelpers import (PlotHelper,\
                           collectiveCollect,\
                           safeCollect,\
                           DDZ,\
                           findLargestRadialGrad)
import numpy as np
from scipy.stats import kurtosis, skew
from scipy.signal import periodogram


# Get radial indices -> equidistanced indices
xind
yind
zind

OR

nXind + center
yind
zind

OR

xind
nYind
zind

# So nYind has higher precedence than yind
# xind etc should still have dimension
# they are called from the outside (e.g. generic driver)
#--------------



# Find the fluctuations in var
self._varAvg   = polAvg(var[tIndSaturatedTurb:, :, :, :])
self._varFluct = var[tIndSaturatedTurb:, :, :, :] - self._varAvg

# Get the units (eventually convert to physical units)
self._var, self.varNormalization, self.varUnits =\
    self.helper.physicalUnitsConverter(var, varName)


# Find the max gradient of the variable (subtracts the guard cells)
_, maxGradInd =\
    findLargestRadialGrad(\
      self._varSteadyState[0:1, :, self.yInd:self.yInd+1, 0:1],\
      dx,\
      self._MXG)

self.timeTraceOfVarFluct\
        ["{},{},{}".format(xInd, actualYInd, zInd)] =\
            self._varFluct[:, xInd, yInd, zInd]
