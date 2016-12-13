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
B0 = (  6.0e-2,   4.0e-2,  2.0e-2)
Lx = (  4.7180,   3.1453,  1.5727)
Ly = (165.1286, 110.0858, 55.0429)

scanParameters  = ("B0", "Lx", "Ly")
series_add = (\
              ("input", "B0", B0),\
              ("geom" , "Lx", Lx),\
              ("geom" , "Ly", Ly),\
             )

directory = "CSDXMagFieldScanAr"

# Set the main options
scanB0.setMainOptions(\
                       directory      = directory       ,\
                       scanParameters = scanParameters  ,\
                       series_add     = series_add      ,\
                       theRunName     = directory       ,\
                       make           = False           ,\
                       boutRunnersNoise = {"vortD":1e-6},\
                     )

# Set the flags
scanB0.setPostProcessingFlags(\
                              justPostProcess            = True ,\
                              postProcessInit            = False,\
                              postProcessExp             = False,\
                              postProcessLin             = False,\
                              postProcessTurb            = True ,\
                              postProcessLinProfiles     = False,\
                              postProcessTurbProfiles    = False,\
                              postProcessProbesAndEnergy = False,\
                              postProcessGrowthRates     = False,\
# FIXME: Look at energy overshoot, and set correct index (starting from linear run) here
                              tIndSaturatedTurb          = None,\
                             )
turbPost=\
  {"driverName" : "single2DDriver" ,\
   "varName"    : "n"              ,\
   "pltName"    : "n"              ,\
   "tSlice"     : slice(0, None, 2),\
   "varyMaxMin" : True             ,\
   "fluctuation": True             ,\
   "mode"       : "perpAndPol"     ,\
  }

scanB0.setTurbulencePostOptions(useDefault             = False,\
                                useFieldPlotterOptions = True ,\
                                useMaxGrad             = True,\
                                **turbPost)

scanB0.setRunOptions(runInit   = True,\
                     runExpand = True,\
                     runLin    = True,\
                     runTurb   = True)

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
