#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

# The options for the run
# =============================================================================
# *****************************************************************************
eiCollisions = [200]
ny = [24]
noise = True
# *****************************************************************************
# Set the temporal domain
restart      = None;
remove_old = False
nout       = [1]
timestep   = [1, 10]
nout *= len(timestep)
directory  = "d-diagnoseNaulin/"
# Shall we make?
make       = False
# =============================================================================


# The options for the post processing function
# =============================================================================
xguards    = False
yguards    = False
xSlice     = 0
ySlice     = 4
zSlice     = 0
showPlot   = False
savePlot   = True
theRunName = "0-d-0-DiagnoseNaulin"
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 48
BOUT_nodes            = 3
BOUT_ppn              = 20
BOUT_walltime         = '00:10:00'
BOUT_run_name         = theRunName
BOUT_queue            = 'xpresq'
post_process_run_name = 'post' + theRunName.capitalize()
# =============================================================================


# Create the runner
# =============================================================================
myRuns = PBS_runner(\
            directory  = directory ,\
            nproc      = nproc ,\
            ny         = ny,\
            # Set temporal domain
            nout       = nout  ,\
            timestep   = timestep,\
            # Copy the source file
            cpy_source = True  ,\
            make       = make  ,\
            restart    = restart,\
            additional = [
                          ('tag',theRunName,0),\
                          ('cst','nuEI',eiCollisions),\
                          ('switch','includeNoise',noise),\
                         ],\
            # PBS options
            BOUT_nodes            = BOUT_nodes           ,\
            BOUT_ppn              = BOUT_ppn             ,\
            BOUT_walltime         = BOUT_walltime        ,\
            BOUT_run_name         = BOUT_run_name        ,\
            BOUT_queue            = BOUT_queue           ,\
            )
# =============================================================================


# Perform the run
# =============================================================================
myRuns.execute_runs(remove_old = remove_old)
# =============================================================================
