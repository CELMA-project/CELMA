#!/usr/bin/env python

"""
Init-file for statistics
"""

from .polAvg import polAvg
from .collectiveCollect import collectiveCollect
from .energy import collectEnergy, plotEnergies
from .probes import PerpPlaneProbes, Probes
from .derivatives import DDZ,\
                         DDX,\
                         findLargestRadialGrad

import matplotlib.cm as cm

colorfunc = cm.viridis
