#!/usr/bin/env python

"""Driver which checks the plots."""

from bout_runners import basic_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotting import combined1D2D

# The options for the run
# =============================================================================
# Set the temporal domain
restart    = None
remove_old = False
nout       = 4
timestep   = 1e-1
internalTime = timestep
directory  = "a-data"
# Shall we make?
make       = True
# The number of processors
nproc = 4
# =============================================================================


# The options for the post processing function
# =============================================================================
xguards       = False
yguards       = False
xSlice        = 0
ySlice        = 8
zSlice        = 0
showPlot      = False
savePlot      = True
useSubProcess = True
theRunName    = "laptop"
# =============================================================================


# Create the runner
# =============================================================================
myRuns = basic_runner(\
                      directory  = directory ,\
                      nproc      = nproc     ,\
                      # Set temporal domain
                      nout       = nout      ,\
                      timestep   = timestep  ,\
                      additional = ('solver', 'timestep', internalTime),\
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
                     post_processing_function = combined1D2D,\
                     # This function will be called every time after
                     # performing a run
                     post_process_after_every_run = True,\
                     # Below are the kwargs arguments being passed to
                     # the post processing function
                     # Switches
                     xguards        = xguards           ,\
                     yguards        = yguards           ,\
                     xSlice         = xSlice            ,\
                     ySlice         = ySlice            ,\
                     zSlice         = zSlice            ,\
                     savePlot       = savePlot          ,\
                     saveFolderFunc = "scanWTagSaveFunc",\
                     theRunName     = theRunName        ,\
                     useSubProcess  = useSubProcess     ,\
                    )
# =============================================================================
