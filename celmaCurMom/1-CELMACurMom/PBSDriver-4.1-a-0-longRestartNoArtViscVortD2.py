#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common/python')
# Sys path is a list of system paths
sys.path.append(commonDir)

from plotting import combinedDriver as postProcess

# The options for the run
# =============================================================================
# *****************************************************************************
includeNoise        = False
forceAddNoise       = False
artViscParVortD     = 0.0
artViscPerpVortD    = 0.0
useHyperViscAzVortD = [True]
# *****************************************************************************
remove_old = False
restart    = "overwrite"
# Uncomment this if you just want to plot
# restart      = None;
restart_from = "a-data/nout_300_timestep_10/nz_128/cst_artViscParVortD_0.0_cst_artViscPerpVortD_0.0_switch_forceAddNoise_False_switch_includeNoise_False_switch_useHyperViscAzVortD_True_tag_3.1-a-0-longRestartNoArtViscVortD_0/"
# Set the spatial domain
nz = 128
# Set the temporal domain
nout       = [500]
timestep   = [30]
directory  = "a-data"
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
tSlice     = slice(-20, None)
showPlot   = False
savePlot   = True
theRunName = "4.1-a-0-longRestartNoArtViscVortD2"
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 24
BOUT_nodes            = 2
BOUT_ppn              = 12
BOUT_walltime         = '48:00:00'
BOUT_run_name         = theRunName
post_process_nproc    = 1
post_process_nodes    = 1
post_process_ppn      = 20
post_process_walltime = '0:29:00'
post_process_queue    = 'xpresq'
post_process_run_name = 'post' + theRunName.capitalize()
# =============================================================================


# Create the runner
# =============================================================================
myRuns = PBS_runner(\
            directory  = directory ,\
            nproc      = nproc ,\
            # Set spatial domain
            nz         = nz,\
            # Set temporal domain
            nout       = nout  ,\
            timestep   = timestep,\
            # Copy the source file
            make       = make  ,\
            restart    = restart,\
            restart_from = restart_from,\
            additional = [
                          ('tag',theRunName,0),\
                          ('switch'      , 'includeNoise'       , includeNoise ),\
                          ('switch'      , 'forceAddNoise'      ,forceAddNoise),\
                          ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                          ('cst'   , 'artViscParVortD'    , artViscParVortD ),\
                          ('cst'   , 'artViscPerpVortD'   , artViscPerpVortD),\
                         ],\
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
                     post_processing_function = postProcess,\
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
                     tSlice         = tSlice            ,\
                     savePlot       = savePlot          ,\
                     saveFolderFunc = "scanWTagSaveFunc",\
                     theRunName     = theRunName        ,\
                    )
# =============================================================================
