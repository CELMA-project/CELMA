#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import GenericScanDriver

# Create object
scanB0 = GenericScanDriver()

# Set the scan
B0 = ( 2.0e-2,)
Lx = ( 1.5727,)
Ly = (55.0429,)

scanParameters  = ("B0", "Lx", "Ly")
series_add = (\
              ("input", "B0", B0),\
              ("geom" , "Lx", Lx),\
              ("geom" , "Ly", Ly),\
             )

directory = "WRONGBCS-CSDXMagFieldScanAr"

# Set the main options
scanB0.setMainOptions(\
                       directory      = directory       ,\
                       scanParameters = scanParameters  ,\
                       series_add     = series_add      ,\
                       theRunName     = directory       ,\
                       make           = False           ,\
                       boutRunnersNoise = {"vort":1e-6} ,\
                     )

# Increase timestep just in case
scanB0.setExpandOptions(BOUT_walltime = "72:00:00")

# Set the flags
scanB0.setPostProcessingFlags(\
                              justPostProcess            = False,\
                              postProcessInit            = False,\
                              postProcessExp             = False,\
                              postProcessLin             = False,\
                              postProcessTurb            = False,\
                              postProcessLinProfiles     = False,\
                              postProcessTurbProfiles    = False,\
                              postProcessProbesAndEnergy = False,\
                              postProcessGrowthRates     = False,\
# FIXME: Look at energy overshoot, and set correct index (starting from linear run) here
                              tIndSaturatedTurb          = None,\
                             )

# Set common plotter options
scanB0.setCommonPlotterOptions(\
                               saveFolderFunc    = "scanWTagSaveFunc",\
                               convertToPhysical = True              ,\
                               showPlot          = False             ,\
                               savePlot          = True              ,\
                               extension         = "png"             ,\
                               useSubProcess     = True              ,\
                              )

# Set probe plotter options
scanB0.setProbePlottersOptions(\
                               nProbes = 5  ,\
                               maxMode = 10 ,\
                               yInd    = 16 ,\
                              )

# Set field plotter options
scanB0.setFieldPlottersOptions(\
                               xguards           = False,\
                               yguards           = False,\
                               xSlice            = 0    ,\
                               ySlice            = 16   ,\
                               zSlice            = 0    ,\
                               axisEqualParallel = False,\
                              )

# Set common runner options
scanB0.setCommonRunnerOptions(\
                              nproc              = 48  ,\
                              cpy_source         = True,\
                              BOUT_nodes         = 3   ,\
                              BOUT_ppn           = 16  ,\
                              post_process_nproc = 1   ,\
                              post_process_nodes = 1   ,\
                              post_process_ppn   = 20  ,\
                             )

# Run
scanB0.runScan()
