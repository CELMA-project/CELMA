#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import GenericScanDriver

# Create object
scanNn = GenericScanDriver()

# Set the scan
nn             = (1.0e15, 5.0e15, 1.0e16, 5.0e16, 1.0e17)
scanParameters = ("nn",)
series_add = ('input', 'nn', nn)

directory = "a1-KiwiFlatNeutralTi0"

# Set the main options
scanNn.setMainOptions(\
                       directory      = directory     ,\
                       scanParameters = scanParameters,\
                       series_add     = series_add    ,\
                       theRunName     = directory     ,\
                       make           = False         ,\
                       varName        = "n"           ,\
                       pltName        = "n"           ,\
                     )

# Set the flags
scanNn.setPostProcessingFlags(\
                              justPostProcess            = True ,\
                              postProcessInit            = True ,\
                              postProcessExp             = True ,\
                              postProcessLin             = True ,\
                              postProcessTurb            = True ,\
                              postProcessLinProfiles     = True ,\
                              postProcessTurbProfiles    = True ,\
                              postProcessProbesAndEnergy = True ,\
                              postProcessGrowthRates     = True ,\
                             )

# Set common plotter options
scanNn.setCommonPlotterOptions(\
                               saveFolderFunc    = "scanWTagSaveFunc",\
                               convertToPhysical = True              ,\
                               showPlot          = False             ,\
                               savePlot          = True              ,\
                               extension         = "png"             ,\
                               useSubProcess     = True              ,\
                              )

# Set probe plotter options
scanNn.setProbePlottersOptions(\
                               nProbes = 5  ,\
                               maxMode = 10 ,\
                               yInd    = 16 ,\
                              )

# Set field plotter options
scanNn.setFieldPlottersOptions(\
                               xguards           = False,\
                               yguards           = False,\
                               xSlice            = 0    ,\
                               ySlice            = 16   ,\
                               zSlice            = 0    ,\
                               axisEqualParallel = False,\
                              )

# Set common runner options
scanNn.setCommonRunnerOptions(\
                              nproc              = 48  ,\
                              cpy_source         = True,\
                              BOUT_nodes         = 3   ,\
                              BOUT_ppn           = 16  ,\
                              post_process_nproc = 1   ,\
                              post_process_nodes = 1   ,\
                              post_process_ppn   = 20  ,\
                             )

# Run
scanNn.runScan()
