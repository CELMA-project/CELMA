#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common/python')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotting import combined1D2D

# The options for the run
# =============================================================================
# Set the temporal domain
restart      = "overwrite"
# Uncomment this if you just want to plot
# restart      = None
restart_from = "../1-CELMA/lessSource/nout_20_timestep_50.0/tag_0-c-0-LessSource_0/"
remove_old = False
timestep   = [5e0]
nout       = [20]*len(timestep)
directory  = "restartBroaderCauchy"
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
saveFolder = "0-d-0-RestartFromCELMA0c0BroaderCauchy"
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 96
BOUT_nodes            = 5
BOUT_ppn              = 20
BOUT_walltime         = '08:00:00'
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
            additional = [('tag',saveFolder,0)           ,\
                         ]                               ,\
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
                     post_processing_function = combined1D2D,\
                     # This function will be called every time after
                     # performing a run
                     post_process_after_every_run = True,\
                     # Below are the kwargs arguments being passed to
                     # the post processing function
                     # Switches
                     xguards    = xguards    ,\
                     yguards    = yguards    ,\
                     xSlice     = xSlice     ,\
                     ySlice     = ySlice     ,\
                     zSlice     = zSlice     ,\
                     savePlot   = savePlot   ,\
                     saveFolder = saveFolder ,\
                    )
# =============================================================================
