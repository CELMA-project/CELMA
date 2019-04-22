#!/usr/bin/env python

"""
Driver which plots the results of the simulations.

**NOTE:** This driver has unfilled fields marked with FIXME
"""

import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from standardPlots import PlotSubmitter

# INPUT
# =============================================================================
# If you would like a mail on finished job enter your mail here
# Example: "john@doe.com"
mail  = None
# If the queuing system uses accounts, set the account here
# Example: "FUA11_SOLF"
account = None
# Usually, the queueing system has its own default queue, if not,
# specify here
# Example: "xfualongprod"
queue = None
# =============================================================================

directory = "BousCSDXMagFieldScanAr"
scanParameter = "B0"

# Create the plotSubmitter
pltSub = PlotSubmitter(directory, scanParameter, boussinesq=True)
pltSub.sub.setMisc(logPath = os.path.join(directory,"postLogs"),\
                   mail    = mail,\
                   account = account)
pltSub.sub.setQueue(queue)

# Set linear slices (found from looking at the Fourier modes)
# FIXME:
tSlices = {\
           "B0_0.02":None,\
           "B0_0.04":None,\
           "B0_0.06":None,\
           "B0_0.08":None,\
           "B0_0.1" :None,\
           }
pltSub.setLinearPhaseTSlices(tSlices)

# Set saturated turbulence slices (found from looking at the energies)
# FIXME:
tSlices = {\
           "B0_0.02":None,\
           "B0_0.04":None,\
           "B0_0.06":None,\
           "B0_0.08":None,\
           "B0_0.1" :None,\
           }
pltSub.setSatTurbTSlices(tSlices)

# Run the post-processing
pltSub.updatePlotSuperKwargs({"extension" : "pdf"})
pltSub.runAnalyticGrowthRates()

pltSub.runBlobs(modes="perp", fluct=True, condition=3)
pltSub.runBlobs(modes="perp", fluct=True, condition=2)
pltSub.runBlobs(modes="perp", fluct=True, condition=4)

pltSub.runBlobDensPDF()

pltSub.runCominedPlots()
pltSub.runEnergy(sliced=False)
pltSub.runEnergy(sliced=True)
pltSub.runFourierModes(sliced=False)
pltSub.runFourierModes(sliced=True)
pltSub.runGrowthRates()
pltSub.runPerformance()
pltSub.runPhaseShift()
pltSub.runPosOfFluct()
pltSub.runPSD2D()
pltSub.runSkewKurt()
pltSub.runSteadyState()
pltSub.runZonalFlow()
pltSub.runTotalFlux()

# Run the animations
pltSub.updatePlotSuperKwargs({"extension" : None})
pltSub.runFields1DAnim()
pltSub.sub.setWalltime("06:00:00")
pltSub.runFields2DAnim(fluct=True)
pltSub.runFields2DAnim(fluct=False)

# Snapshots plot
# Obtain frames to see evolution
# FIXME: (Maybe ok values to use, but check this)
# NOTE: Plot 4,6 and 8 are chosen
# B=0.01
start   = 2055
end     = 2155
pics    = 10
frameNr = np.linspace(start, end, pics)
# Run runSnapShotsSameScanVal without vMaxVmin in order to see maxMin
# FIXME:
maxMin = (1.35e18, 7.50e18)
maxMin = None
turbSlices = tuple(slice(int(frame),int(frame)) for frame in frameNr)
pltSub.runSnapShotsSameScanVal("param0",turbSlices,fluct=False,vMaxVMin=maxMin)
