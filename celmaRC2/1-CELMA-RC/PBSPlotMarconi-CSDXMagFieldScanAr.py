#!/usr/bin/env python

"""Driver which plots the results of the simulations."""

import numpy as np

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
pltSub.sub.setMisc(logPath = os.path.join(directory,"postLogs"),\
                   mail    = "mmag@fysik.dtu.dk",\
                   account = "FUA11_SOLF")
pltSub.sub.setQueue("xfualongprod")

# Set linear slices
tSlices = {\
           "B0_0.02":slice(80,300)  ,\
           "B0_0.04":slice(800,1250),\
           "B0_0.06":slice(180,280) ,\
           "B0_0.08":slice(100,225) ,\
           "B0_0.1" :slice(80,210)  ,\
          }
pltSub.setLinearPhaseTSlices(tSlices)

# Set saturated turbulence slices
tSlices = {\
           "B0_0.02":None            ,\
           "B0_0.04":None            ,\
           "B0_0.06":slice(1500,None),\
           "B0_0.08":slice(1200,None),\
           "B0_0.1" :slice(1000,None),\
           }
pltSub.setSatTurbTSlices(tSlices)

# Run the post-processing
pltSub.updatePlotSuperKwargs({"extension" : "pdf"})
pltSub.runAnalyticGrowthRates()
pltSub.runCominedPlots()
pltSub.runEnergy(sliced=False)
pltSub.runEnergy(sliced=True)
pltSub.runFourierModes(sliced=False)
pltSub.runFourierModes(sliced=True)
pltSub.runGrowthRates()
pltSub.runMagnitudeSpectrum()
pltSub.runPerformance(allFolders=False)
pltSub.runPerformance(allFolders=True)
pltSub.runPhaseShift()
pltSub.runPosOfFluct()
pltSub.runPSD2D()
pltSub.runSkewKurt()
pltSub.runSteadyState()
pltSub.runZonalFlow()
# Post processing taking longer time
pltSub.sub.setWalltime("00:30:00")
pltSub.runTotalFlux()

# Run the animations
pltSub.updatePlotSuperKwargs({"extension" : None})
pltSub.sub.setWalltime("01:00:00")
pltSub.runFields1DAnim()
pltSub.sub.setWalltime("06:00:00")
pltSub.runFields2DAnim(fluct=True)
pltSub.runFields2DAnim(fluct=False)

# Snapshots plot
# Obtain evolution of the mode
modeSlices = (\
              slice(168, 168),\
              slice(183, 183),\
              slice(198, 198),\
              )
pltSub.runSnapShotsSameScanVal("param0", modeSlices, fluct=True, yInd=50)

# Obtain the different modes
modesSlices = {\
               "B0_0.02":slice( 147,  147),\
               "B0_0.04":slice(1500, 1500),\
               "B0_0.06":slice( 300,  300),\
               "B0_0.08":slice( 230,  230),\
               "B0_0.1" :slice( 198,  198),\
               }
pltSub.runSnapShotDifferentScanVals(modesSlices, fluct=True, yInd=50)

# Obtain frames to see evolution
# NOTE: Plot 4,6 and 8 are chosen
# B=0.01
start   = 2055
end     = 2155
pics    = 10
frameNr = np.linspace(start, end, pics)
# Run runSnapShotsSameScanVal without vMaxVmin in order to see maxMin
maxMin = (1.35e18, 7.50e18)
turbSlices = tuple(slice(int(frame),int(frame)) for frame in frameNr)
pltSub.runSnapShotsSameScanVal("param0",turbSlices,fluct=False,vMaxVMin=maxMin)

# Obtain the turbulence fluctuations
# FIXME:
turbSlices = (\
              slice(2122, 2122),\
             )
pltSub.runSnapShotsSameScanVal("param0", turbSlices, fluct=True)
