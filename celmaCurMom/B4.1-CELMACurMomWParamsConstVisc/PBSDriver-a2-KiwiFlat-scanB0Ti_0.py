#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import GenericScanDriver

# Create object
scanB0 = GenericScanDriver()

# Set the scan
B0 = (1.0e-1  , 9.0e-2  , 8.0e-2  , 7.0e-2  , 6.0e-2  , 5.0e-2  )
Lx = (6.5393  , 5.8854  , 5.2315  , 4.5775  , 3.9236  , 3.2697  )
Ly = (366.2017, 329.5815, 292.9613, 256.3412, 219.7210, 183.1008)
scanParameters  = ("B0", "Lx", "Ly")
series_add = (\
              ('input', 'B0', B0),\
              ('geom' , 'Lx', Lx),\
              ('geom' , 'Ly', Ly),\
             )

directory = "a2-KiwiFlatMagFieldTi0"

# Set the main options
scanB0.setMainOptions(\
                       directory      = directory     ,\
                       scanParameters = scanParameters,\
                       series_add     = series_add    ,\
                       theRunName     = directory     ,\
                       make           = False         ,\
                       varName        = "n"           ,\
                       pltName        = "n"           ,\
                     )

# Set the flags
scanB0.setPostProcessingFlags(\
                              justPostProcess            = True ,\
                              postProcessInit            = False,\
                              postProcessExp             = False,\
                              postProcessLin             = False,\
                              postProcessTurb            = False,\
                              postProcessLinProfiles     = False,\
                              postProcessTurbProfiles    = False,\
                              postProcessProbesAndEnergy = True ,\
                              postProcessGrowthRates     = False,\
                              # Calculated from the energy overshoot
                              # FIXME: Check when this happens
                              tIndSaturatedTurb          = None ,\
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
