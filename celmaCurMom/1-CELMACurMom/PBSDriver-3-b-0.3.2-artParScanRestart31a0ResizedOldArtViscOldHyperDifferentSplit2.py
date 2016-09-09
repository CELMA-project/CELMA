#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotting import combinedDriver as postProcess

# The options for the run
# =============================================================================
# *****************************************************************************
includeNoise        = False
forceAddNoise       = False
artPerp             = [2.0e-1]
artPar              = [10.0, 4.0, 4.0e-1]
useHyperViscAzVortD = [True]
artHyperAzVortD     = 16
# *****************************************************************************
remove_old = False
restart    = "overwrite"
# Uncomment this if you just want to plot
# restart      = None;
restart_from = "b-noArtVortD/nout_0_timestep_1e-10/nx_34_ny_48_nz_256/cst_artPar_2.5_switch_forceAddNoise_False_switch_includeNoise_False_switch_useHyperViscAzVortD_True_tag_3-b-0.0.3.2-31a0Resized_0/"
# Set the spatial domain
nx = 16*2 + 2
ny = 24*2
nz = 128*2
# Set the temporal domain
nout       = [101]
timestep   = [30]
directory  = "b-noArtVortD"
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
tSlice     = slice(-20, None)
showPlot   = False
savePlot   = True
theRunName = "3-b-0.3.2-artParScanRestart31a0ResizedOldArtViscOldHyperDifferentSplit2"
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 64
BOUT_nodes            = 4
BOUT_ppn              = 16
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
            nx         = nx,\
            ny         = ny,\
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
                          ('switch', 'includeNoise'       , includeNoise ),\
                          ('switch', 'forceAddNoise'      ,forceAddNoise),\
                          ('switch', 'useHyperViscAzVortD',useHyperViscAzVortD),\
                          ('cst'   , 'artPerp'            ,artPerp),\
                          ('cst'   , 'artPar'             ,artPar),\
                          ('cst'   , 'arthyperazvortd' ,artHyperAzVortD),\
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
