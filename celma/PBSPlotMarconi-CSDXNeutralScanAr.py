#!/usr/bin/env python

"""Driver which plots the results of the simulations."""

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

directory = "CSDXNeutralScanAr"
scanParameter = "nn"

# Create the plotSubmitter
pltSub = PlotSubmitter(directory, scanParameter)
pltSub.sub.setMisc(logPath = os.path.join(directory,"postLogs"),\
                   mail    = mail,\
                   account = account)
pltSub.sub.setQueue(queue)

# Set linear slices (found from looking at the Fourier modes)
tSlices = {\
    "nn_2.5e+18"               : slice(50, 300),\
    "nn_6.666666666666668e+18" : slice(50, 225),\
    "nn_1.5e+19"               : slice(50, 210),\
    "nn_4e+19"                 : slice(50, 210),\
    "nn_9.9e+20"               : slice(50, 500),\
          }
pltSub.setLinearPhaseTSlices(tSlices)

# Set saturated turbulence slices (found from looking at the energies)
tSlices = {\
    "nn_2.5e+18"               : slice(int(1.5e3), None),\
    "nn_6.666666666666668e+18" : slice(int(1.5e3), None),\
    "nn_1.5e+19"               : slice(int(1.5e3), None),\
    "nn_4e+19"                 : slice(int(1.5e3), None),\
    "nn_9.9e+20"               : slice(int(4.5e3), None),\
           }
pltSub.setSatTurbTSlices(tSlices)

# Run the post-processing
pltSub.updatePlotSuperKwargs({"extension" : "pdf"})

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
pltSub.runKThetaSpectrum()
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
# Obtain evolution of the mode
# FIXME
modeSlices = (\
              slice(0, 0),\
              slice(0, 0),\
              slice(0, 0),\
              )
pltSub.runSnapShotsSameScanVal("param0", modeSlices, fluct=True, yInd=50)

# Obtain the different modes
# FIXME
modesSlices = {\
               "nn_2.5e+18"               : None,\
               "nn_6.666666666666668e+18" : None,\
               "nn_1.5e+19"               : None,\
               "nn_4e+19"                 : None,\
               "nn_9.9e+20"               : None,\
               }
pltSub.runSnapShotDifferentScanVals(modesSlices, fluct=True, yInd=50)

# Obtain frames to see evolution
# FIXME
# NOTE: Plot 5,6 and 7 are chosen
# nn_9.9e+20
start   = 2055
end     = 2155
pics    = 10
frameNr = np.linspace(start, end, pics)
# Run runSnapShotsSameScanVal without vMaxVmin in order to see maxMin
# FIXME
maxMin = (1.35e18, 7.50e18)
turbSlices = tuple(slice(int(frame),int(frame)) for frame in frameNr)
pltSub.runSnapShotsSameScanVal("param0",turbSlices,fluct=False,vMaxVMin=maxMin)

# Obtain the turbulence fluctuations
# FIXME
turbSlices = (\
              slice(2122, 2122),\
             )
pltSub.runSnapShotsSameScanVal("param0", turbSlices, fluct=True)
