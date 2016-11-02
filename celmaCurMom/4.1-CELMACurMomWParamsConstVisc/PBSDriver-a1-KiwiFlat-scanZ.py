#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import GenericScanDriver

# Create object
scanZ = GenericScanDriver()

# Set the scan
length = (1       , 2       , 4       , 6       , 8       , 10       )
Ly     = (102.2235, 204.4469, 408.8938, 613.3408, 817.7877, 1022.2346)
scanParameters  = ("length", "Ly")
series_add = (\
              ('input', 'length', length),\
              ('geom' , 'Ly'    , Ly),\
             )

directory = "a1-KiwiFlatZ"

# Set the main options
scanZ.setMainOptions(\
                       directory      = directory     ,\
                       scanParameters = scanParameters,\
                       series_add     = series_add    ,\
                       theRunName     = directory     ,\
                       make           = False         ,\
                       varName        = "n"           ,\
                       pltName        = "n"           ,\
                     )

# Set the flags
scanZ.setPostProcessingFlags(\
                              justPostProcess            = True ,\
                              postProcessInit            = False,\
                              postProcessExp             = False,\
                              postProcessLin             = False,\
                              postProcessTurb            = False,\
                              postProcessLinProfiles     = False,\
                              postProcessTurbProfiles    = False,\
                              postProcessProbesAndEnergy = True ,\
                              postProcessGrowthRates     = True ,\
                              # FIXME: Check that this is true
                              # Calculated from the energy overshoot
                              tIndSaturatedTurb          = None ,\
                             )

# Set common plotter options
scanZ.setCommonPlotterOptions(\
                               saveFolderFunc    = "scanWTagSaveFunc",\
                               convertToPhysical = True              ,\
                               showPlot          = False             ,\
                               savePlot          = True              ,\
                               extension         = "png"             ,\
                               useSubProcess     = True              ,\
                              )

# Set probe plotter options
scanZ.setProbePlottersOptions(\
                               nProbes = 5  ,\
                               maxMode = 10 ,\
                               yInd    = 16 ,\
                              )

# Set field plotter options
scanZ.setFieldPlottersOptions(\
                               xguards           = False,\
                               yguards           = False,\
                               xSlice            = 0    ,\
                               ySlice            = 16   ,\
                               zSlice            = 0    ,\
                               axisEqualParallel = False,\
                              )

# Set common runner options
scanZ.setCommonRunnerOptions(\
                              nproc              = 48  ,\
                              cpy_source         = True,\
                              BOUT_nodes         = 3   ,\
                              BOUT_ppn           = 16  ,\
                              post_process_nproc = 1   ,\
                              post_process_nodes = 1   ,\
                              post_process_ppn   = 20  ,\
                             )

# Run
scanZ.runScan()
