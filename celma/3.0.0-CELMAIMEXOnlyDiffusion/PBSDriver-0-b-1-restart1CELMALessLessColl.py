#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
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
eiCollisions = [550, 500]
# *****************************************************************************
# Set the temporal domain
restart      = "overwrite"
# Uncomment this if you just want to plot
# restart      = None;
restart_from = "/mnt/mimer/mmag/CELMA-dev/celma/1-CELMA/lessSource/nout_20_timestep_50.0/tag_0-c-0-LessSource_0"
remove_old   = False
nout       = [20]
timestep   = [5e0]
internalTime = [5e-4]
directory  = "b-restartFrom1CELMA"
# Shall we make?
make       = False
# =============================================================================


# The options for the post processing function
# =============================================================================
xguards    = False
yguards    = False
xSlice     = 0
ySlice     = 8
zSlice     = 0
showPlot   = False
savePlot   = True
saveFolder = "0-b-1-restart1CELMALessLessColl"
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 96
BOUT_nodes            = 5
BOUT_ppn              = 20
BOUT_walltime         = '06:00:00'
BOUT_run_name         = saveFolder
post_process_nproc    = 1
post_process_nodes    = 1
post_process_ppn      = 20
post_process_walltime = '0:29:00'
post_process_queue    = 'xpresq'
post_process_run_name = 'post' + saveFolder.capitalize()
# =============================================================================


# Create the runner
# =============================================================================
myRuns = PBS_runner(\
            directory  = directory ,\
            nproc      = nproc ,\
            # Set temporal domain
            nout       = nout  ,\
            timestep   = timestep,\
            # Copy the source file
            cpy_source = True  ,\
            make       = make  ,\
            restart    = restart,\
            restart_from = restart_from,\
            # Tag (used to catalogize the runs)
            additional = [\
                          ('tag',saveFolder,0)                ,\
                          ('solver', 'timestep', internalTime),\
                          ('cst','nuEI',eiCollisions),\
                         ]                      ,\
            # PBS options
            BOUT_nodes            = BOUT_nodes           ,\
            BOUT_ppn              = BOUT_ppn             ,\
            BOUT_walltime         = BOUT_walltime        ,\
            BOUT_run_name         = BOUT_run_name        ,\
            post_process_nproc    = post_process_nproc   ,\
            post_process_nodes    = post_process_nodes   ,\
            post_process_ppn      = post_process_ppn     ,\
            post_process_walltime = post_process_walltime,\
            post_process_queue    = post_process_queue   ,\
            post_process_run_name = post_process_run_name,\
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
                     xguards    = xguards    ,\
                     yguards    = yguards    ,\
                     xSlice     = xSlice     ,\
                     ySlice     = ySlice     ,\
                     zSlice     = zSlice     ,\
                     savePlot   = savePlot   ,\
                     saveFolder = saveFolder ,\
                    )
# =============================================================================
