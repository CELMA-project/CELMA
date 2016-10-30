#!/usr/bin/env python

"""Driver which checks the plots."""

from bout_runners import basic_runner

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import postBoutRunner

# The options for the run
# =============================================================================
# Set the temporal domain
restart    = None
remove_old = False
nout       = (2,)
timestep   = (1e-10,)
directory  = "a1-KiwiFlatMagField"
# Shall we make?
make       = False
# The number of processors
nproc = 4
# =============================================================================


# The options for the post processing function
# =============================================================================
xguards           = False
yguards           = False
xSlice            = 0
ySlice            = 8
zSlice            = 0
tSlice            = None
showPlot          = False
savePlot          = True
theRunName        = "test"
saveFolderFunc    = "scanWTagSaveFunc"
convertToPhysical = False
extension         = "png"
axisEqualParallel = False
useSubProcess     = True
# =============================================================================


# Create the runner
# =============================================================================
myRuns = basic_runner(\
                      directory  = directory ,\
                      nproc      = nproc     ,\
                      # Set temporal domain
                      nout       = nout      ,\
                      timestep   = timestep  ,\
                      # Copy the source file
                      cpy_source = True      ,\
                      make       = make      ,\
                      restart    = restart   ,\
                      )
# =============================================================================


# Perform the run
# =============================================================================
myRuns.execute_runs(\
                     remove_old               = remove_old    ,\
                     post_processing_function = postBoutRunner,\
                     # This function will be called every time after
                     # performing a run
                     post_process_after_every_run = True,\
                     # Below are the kwargs arguments being passed to
                     # the post processing function
                     driverName     = "plot1D2DAndFluctDriver",\
                     tSlice         = tSlice             ,\
                     theRunName     = theRunName         ,\
                     # Extra kwargs
                     saveFolderFunc   = saveFolderFunc   ,\
                     convertToPhysical= convertToPhysical,\
                     showPlot         = showPlot         ,\
                     savePlot         = savePlot         ,\
                     extension        = extension        ,\
                     useSubProcess    = useSubProcess    ,\
                     xguards          = xguards          ,\
                     yguards          = yguards          ,\
                     xSlice           = xSlice           ,\
                     ySlice           = ySlice           ,\
                     zSlice           = zSlice           ,\
                     axisEqualParallel= axisEqualParallel,\
                    )
# =============================================================================
