#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import basic_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import postBoutRunner

# The options for the run
# =============================================================================
# *****************************************************************************
ownFilterType = "none"
# *****************************************************************************
# Set the spatial domain
nz = 16
nx = 32 + 2
ny = 48
nproc = 4
# Set the temporal domain
restart    = "overwrite"
restart_from = "test-resize/nout_0_timestep_1e-09/nz_16/ownFilters_type_none_tag_driverResizeTest3_0/"
remove_old = False
timestep   = [1e-9]
nout       = [0]
directory  = "test-resize"
# Shall we make?
make       = False
# =============================================================================


# The options for the post processing function
# =============================================================================
xguards    = False
yguards    = False
xSlice     = 0
ySlice     = 8*2
zSlice     = 0
showPlot   = False
savePlot   = True
theRunName = "driverResizeTest7"
# =============================================================================


# Create the runner
# =============================================================================
myRuns = basic_runner(\
            restart_from = restart_from,\
            directory  = directory ,\
            nproc      = nproc ,\
            # Set spatial domain
            nx         = nx,\
            nz         = nz,\
            ny         = ny,\
            # Set temporal domain
            nout       = nout  ,\
            timestep   = timestep,\
            # Copy the source file
            make       = make  ,\
            restart    = restart,\
            additional = [
                          ('tag',theRunName,0),\
                          ('ownFilters'  , 'type', ownFilterType),\
                         ]
                      )
# =============================================================================


# Perform the run
# =============================================================================
myRuns.execute_runs(\
                     remove_old               = remove_old,\
                     post_processing_function = postBoutRunner,\
                     # This function will be called every time after
                     # performing a run
                     post_process_after_every_run = True,\
                     # Below are the kwargs arguments being passed to
                     # the post processing function
                     # Switches
                     driverName        = "plot1D2DAndFluctDriver",\
                     xguards        = xguards           ,\
                     yguards        = yguards           ,\
                     xSlice         = xSlice            ,\
                     ySlice         = ySlice            ,\
                     zSlice         = zSlice            ,\
                     savePlot       = savePlot          ,\
                     saveFolderFunc = "scanWTagSaveFunc",\
                     theRunName     = theRunName        ,\
                    )
# =============================================================================
