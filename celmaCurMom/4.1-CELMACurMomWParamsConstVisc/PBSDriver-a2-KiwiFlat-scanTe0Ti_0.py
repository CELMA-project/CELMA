#!/usr/bin/env python

"""Driver which runs using PBS."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import GenericScanDriver

# Create object
scanTe0 = GenericScanDriver()

# Set the scan
Te0 = (14.0    , 12.0    , 10.0    , 8.0     , 6.0     , 5.0     )
Lx  = (4.8296  , 5.0295  , 5.2565  , 5.5174  , 5.8213  , 5.9934  )
Ly  = (270.4579, 281.6539, 294.3664, 308.9719, 325.9921, 335.6294)
scanParameters  = ("Te0", "Lx", "Ly")
series_add = (\
              ('input', 'Te0', Te0),\
              ('geom' , 'Lx' , Lx),\
              ('geom' , 'Ly' , Ly),\
             )

directory = "a2-KiwiFlatElTempTi0"

# Set the main options
scanTe0.setMainOptions(\
                       directory      = directory     ,\
                       scanParameters = scanParameters,\
                       series_add     = series_add    ,\
                       theRunName     = directory     ,\
                       make           = False         ,\
                       varName        = "n"           ,\
                       pltName        = "n"           ,\
                     )

# Set the flags
scanTe0.setPostProcessingFlags(\
                              justPostProcess            = False,\
                              postProcessInit            = False,\
                              postProcessExp             = False,\
                              postProcessLin             = False,\
                              postProcessTurb            = False,\
                              postProcessLinProfiles     = False,\
                              postProcessTurbProfiles    = False,\
                              postProcessProbesAndEnergy = False,\
                              postProcessGrowthRates     = False,\
                              # FIXME: Check that this is true
                              # Calculated from the energy overshoot
                              tIndSaturatedTurb          = 600  ,\
                             )

# Set common plotter options
scanTe0.setCommonPlotterOptions(\
                               saveFolderFunc    = "scanWTagSaveFunc",\
                               convertToPhysical = True              ,\
                               showPlot          = False             ,\
                               savePlot          = True              ,\
                               extension         = "png"             ,\
                               useSubProcess     = True              ,\
                              )

# Set probe plotter options
scanTe0.setProbePlottersOptions(\
                               nProbes = 5  ,\
                               maxMode = 10 ,\
                               yInd    = 16 ,\
                              )

# Set field plotter options
scanTe0.setFieldPlottersOptions(\
                               xguards           = False,\
                               yguards           = False,\
                               xSlice            = 0    ,\
                               ySlice            = 16   ,\
                               zSlice            = 0    ,\
                               axisEqualParallel = False,\
                              )

# Set common runner options
scanTe0.setCommonRunnerOptions(\
                              nproc              = 48  ,\
                              cpy_source         = True,\
                              BOUT_nodes         = 3   ,\
                              BOUT_ppn           = 16  ,\
                              post_process_nproc = 1   ,\
                              post_process_nodes = 1   ,\
                              post_process_ppn   = 20  ,\
                             )

# Run
scanTe0.runScan()
