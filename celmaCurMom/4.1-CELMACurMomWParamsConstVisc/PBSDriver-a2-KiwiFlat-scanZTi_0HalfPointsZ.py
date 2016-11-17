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
length = (4       , 5       , 6       , 8        , 10       ,  12       )
Ly     = (553.6449, 692.0561, 830.4673, 1107.2898, 1384.1122,  1660.9347)
scanParameters  = ("len", "Ly")
series_add = (\
              ('input', 'length', length),\
              ('geom' , 'Ly'    , Ly),\
             )

directory = "a2-KiwiFlatZTi0HalfPointsZ"

# Set the main options
scanZ.setMainOptions(\
                       directory             = directory     ,\
                       scanParameters        = scanParameters,\
                       series_add            = series_add    ,\
                       theRunName            = directory     ,\
                       make                  = False         ,\
                       timeStepMultiplicator = 10            ,\
                     )

# Set the flags
scanZ.setPostProcessingFlags(\
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
                              tIndSaturatedTurb          = None ,\
                             )

# Modify nz from default
scanZ.setExpandOptions(nz = 128)

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
