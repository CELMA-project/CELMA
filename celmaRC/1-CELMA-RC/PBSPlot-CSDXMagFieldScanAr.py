#!/usr/bin/env python

"""Driver which plots the results of the simulations."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from standardPlots import PlotSubmitter

directory = "CSDXMagFieldScanAr"
scanParameter = "B0"

# Create the plotSubmitter
pltSub = PlotSubmitter(directory, scanParameter)

# Run the fourier modes
# pltSub.runFourierModes(sliced=False)
# Set linear slices and plot the sliced fourier modes
tSlices = {\
           "B0_0.02":slice(80,240)  ,\
           "B0_0.04":slice(800,1250),\
           "B0_0.06":slice(180,300) ,\
           "B0_0.08":slice(100,225) ,\
           "B0_0.1" :slice(80,210)  ,\
           }
plt.setLinearPhaseTSlices(tSlices)

# pltSub.runFourierModes(sliced=True)
# pltSub.runGrowthRates()
# pltSub.runEnergy(sliced=False)
tSlices = {\
           "B0_0.02":None,\
           "B0_0.04":None,\
           "B0_0.06":slice(1200,None),\
           "B0_0.08":slice(1000,None),\
           "B0_0.1" :None,\
           }
plt.setSaturatedTurbulenceTSlices(tSlices)
# pltSub.runEnergy(sliced=True)
# pltSub.runPosOfFluct()
# pltSub.runPSD2D()
# pltSub.runSkewKurt()
# pltSub.runZonalFlow()
# pltSub.runCominedPlots()
pltSub.runPerformance(allFolders=False)
pltSub.runPerformance(allFolders=True)
