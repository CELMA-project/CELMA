#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common/python')
# Sys path is a list of system paths
sys.path.append(commonDir)

from postProcessing.plotting import combinedDriver

# The options for the run
# =============================================================================
# *****************************************************************************
artPar = [5e1, 5e0, 5e-1, 5e-2, 5e-3, 5e-4]
# *****************************************************************************
remove_old = False
restart    = "overwrite"
# Uncomment this if you just want to plot
# restart      = None;
restart_from = "../8-CELMASplit/f-moreSource3Brackets/nout_300_timestep_10/nz_128/cst_artHyperAzVortD_1.0_ownFilters_type_radialLowPass_ownOperators_type_BasicBrackets_switch_forceAddNoise_True_switch_includeNoise_True_switch_saveDdt_True_switch_useHyperViscAzVortD_True_tag_2-e-3.1-moreSAddNRadialLPHyperVBBrackets_0_theSource_a_0.02/"
# Set the spatial domain
nz = 128
# Set the temporal domain
nout       = [20]
timestep   = [1]
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
tSlice     = None
showPlot   = False
savePlot   = True
theRunName = "a-1-R8-2e31-parViscScan"
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 24
BOUT_nodes            = 2
BOUT_ppn              = 12
BOUT_walltime         = '03:00:00'
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
                          ('cst','artPar',artPar),\
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
                     post_processing_function = combinedDriver,\
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
