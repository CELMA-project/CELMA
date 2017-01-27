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
pltSub.runPosOfFluct()
pltSub.runPSD2D()
pltSub.runSkewKurt()
pltSub.runSteadyState()
pltSub.runZonalFlow()
# Post processing taking longer time
# FIXME:
pltSub.sub.setQueue("workq")
pltSub.sub.setWalltime("00:30:00")
pltSub.runTotalFlux()

# Run the animations
pltSub.updatePlotSuperKwargs({"extension" : None})
pltSub.sub.setWalltime("01:00:00")
pltSub.runFields1DAnim()
# FIXME:
pltSub.sub.setQueue("fatq")
pltSub.sub.setWalltime("06:00:00")
pltSub.runFields2DAnim(fluct=True)
pltSub.runFields2DAnim(fluct=False)

# Snapshots plot
# Obtain evolution of the mode
# FIXME:
modeSlices = (\
              slice(None, None),\
              slice(None, None),\
              slice(None, None),\
             )
pltSub.runSnapShotsSameScanVal("param0", modeSlices, fluct=True)

# Obtain the different modes
# FIXME:
modesSlices = {\
               "B0_0.02":slice(None, None),\
               "B0_0.04":slice(None, None),\
               "B0_0.06":slice(None, None),\
               "B0_0.08":slice(None, None),\
               "B0_0.1" :slice(None, None),\
               }
pltSub.runSnapShotDifferentScanVals(modesSlices, fluct=True)

# Obtain frames to see evolution
# FIXME:
# B=0.08
start   = None
end     = None
pics    = 3
frameNr = np.linspace(start, end, pics)
# Run runSnapShotsSameScanVal without vMaxVmin in order to see maxMin
# FIXME:
maxMin = (None, None)
turbSlices = tuple(slice(int(frame),int(frame)) for frame in frameNr)
pltSub.runSnapShotsSameScanVal("param1",turbSlices,fluct=False,vMaxVMin=maxMin)

# Obtain the turbulence fluctuations
# FIXME:
turbSlices = (\
              slice(None, None),\
             )
pltSub.runSnapShotsSameScanVal("param0", turbSlices, fluct=True)
