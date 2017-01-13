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
pltSub.sub.setMisc(logPath = os.path.join(directory,"postLogs"))

# Set linear slices
tSlices = {\
           "B0_0.02":slice(80,240)  ,\
           "B0_0.04":slice(800,1250),\
           "B0_0.06":slice(180,300) ,\
           "B0_0.08":slice(100,225) ,\
           "B0_0.1" :slice(80,210)  ,\
           }
pltSub.setLinearPhaseTSlices(tSlices)

# Set saturated turbulence slices
tSlices = {\
           "B0_0.02":None,\
           "B0_0.04":None,\
           "B0_0.06":slice(1200,None),\
           "B0_0.08":slice(1000,None),\
           "B0_0.1" :None,\
           }
pltSub.setSatTurbTSlices(tSlices)

# Run the post-processing
pltSub.updatePlotSuperKwargs({"extension" : "pdf"})
pltSub.runCominedPlots()
pltSub.runEnergy(sliced=False)
pltSub.runEnergy(sliced=True)
pltSub.runFourierModes(sliced=False)
pltSub.runFourierModes(sliced=True)
pltSub.runGrowthRates()
pltSub.runPerformance(allFolders=False)
pltSub.runPerformance(allFolders=True)
pltSub.runPosOfFluct()
pltSub.runPSD2D()
pltSub.runSkewKurt()
pltSub.runSteadyState()
pltSub.runZonalFlow()

# Run the animations
pltSub.updatePlotSuperKwargs({"extension" : None})
pltSub.sub.setWalltime("01:00:00")
pltSub.runFields1DAnim()
pltSub.sub.setQueue("fatq")
pltSub.sub.setWalltime("06:00:00")
pltSub.runFields2DAnim(fluct=True)
pltSub.runFields2DAnim(fluct=False)

# Snapshots plot
# Obtain evolution of the mode
modeSlices = (\
              slice(138, 138),\
              slice(162, 162),\
              slice(184, 184),\
             )
pltSub.runSnapShotsSameScanVal("param0", modeSlices, fluct=True)

# Obtain the different modes
modesSlices = {\
               "B0_0.02":slice( 147,  147),\
               "B0_0.04":slice(1500, 1500),\
               "B0_0.06":slice( 300,  300),\
               "B0_0.08":slice( 230,  230),\
               "B0_0.1" :slice( 152,  152),\
               }
pltSub.runSnapShotDifferentScanVals(modesSlices, fluct=True)

# Obtain the turbulence
turbSlices = (\
              # Steady state
              slice(  0,   0),\
              # Violent overshoot
              slice(464, 464),\
              # Plasma off center
              slice(862, 862),\
             )
pltSub.runSnapShotsSameScanVal("param0", turbSlices, fluct=False)

# Obtain the turbulence fluctuations
turbSlices = (\
              slice(961, 961),\
             )
pltSub.runSnapShotsSameScanVal("param0", turbSlices, fluct=True)
