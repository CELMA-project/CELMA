#!/usr/bin/env python

"""Driver which checks the plots."""

from bout_runners import basic_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common/python')
# Sys path is a list of system paths
sys.path.append(commonDir)

from plotting import parPerpDriver

# The options for the run
# =============================================================================
width = [500]
width = [el*2 for el in width]
# Set the temporal domain
restart    = None
remove_old = True
nout       = 2
timestep   = 1e-10
directory  = "fullSOriginalParams"
# Shall we make?
make       = False
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
saveFolder    = "driverTest"
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
                      additional = [
                                    ('theSource','width',width),\
                                   ],\
                      )
# =============================================================================


# Perform the run
# =============================================================================
myRuns.execute_runs(\
                     remove_old               = remove_old    ,\
                     post_processing_function = parPerpDriver     ,\
                     # This function will be called every time after
                     # performing a run
                     post_process_after_every_run = True,\
                     # Below are the kwargs arguments being passed to
                     # the post processing function
                     # Switches
                     xguards       = xguards      ,\
                     yguards       = yguards      ,\
                     xSlice        = xSlice       ,\
                     ySlice        = ySlice       ,\
                     zSlice        = zSlice       ,\
                     savePlot      = savePlot     ,\
                     showPlot      = showPlot     ,\
                     saveFolder    = saveFolder   ,\
                     useSubProcess = useSubProcess,\
                    )
# =============================================================================
