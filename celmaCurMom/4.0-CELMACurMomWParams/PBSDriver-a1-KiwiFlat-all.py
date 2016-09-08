#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common/python')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.plotting import combined1D2D as postProcess
from CELMAPython.plotting import combinedDriver as postProcess2

# If you just want to post-process
justPostProcess = False

# Common options
# =============================================================================
# The options for the post processing function
# *****************************************************************************
xguards    = False
yguards    = False
xSlice     = 0
ySlice     = 8*2
zSlice     = 0
showPlot   = False
savePlot   = True
# *****************************************************************************

# Constructor options
# *****************************************************************************
remove_old = False
directory  = "newKiwiFlat"
make       = False
# *****************************************************************************

# The PBS options
# *****************************************************************************
# Specify the numbers used for the BOUT runs
nproc                 = 48
BOUT_nodes            = 3
BOUT_ppn              = 16
post_process_nproc    = 1
post_process_nodes    = 1
post_process_ppn      = 20
post_process_walltime = '0:29:00'
post_process_queue    = 'xpresq'
# *****************************************************************************
# =============================================================================


# Init and expand options
# =============================================================================
# Filter
ownFilterType = "none"
# Spatial domain
nout       = [2]
# Post processing option
tSlice     = slice(-2, None)
# =============================================================================


# Init runner
# =============================================================================
# Init options
# *****************************************************************************
# Name
theRunName = "a1-KiwiFlat-0-initialize"
# Set the spatial domain
nz = 1
# Set the temporal domain
restart    = None
timestep   = [2e3]
# Specify the numbers used for the BOUT runs
BOUT_walltime         = '03:00:00'
BOUT_run_name         = theRunName
post_process_run_name = 'post' + theRunName.capitalize()
# *****************************************************************************

initRunner = PBS_runner(\
                directory  = directory ,\
                nproc      = nproc ,\
                # Set spatial domain
                nz         = nz,\
                # Set temporal domain
                nout       = nout  ,\
                timestep   = timestep,\
                # Copy the source file
                cpy_source = True  ,\
                make       = make  ,\
                restart    = restart,\
                additional = [
                              ('tag',theRunName,0),\
                              ('ownFilters'  , 'type', ownFilterType),\
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

dmp_folders, PBS_ids = initRunner.execute_runs(\
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

if justPostProcess:
    restart = None

# Expand, Linear and turb options
# =============================================================================
# Set the temporal domain
restart    = "overwrite"
# =============================================================================


# The expand runner
# =============================================================================
# Expand runner options
# *****************************************************************************
# Set the spatial domain
nz = 256
# Set the temporal domain
timestep   = [50]
# From previous outputs
restart_from = dmp_folders[0]
# Name
theRunName = "a1-KiwiFlat-1-expand"
# PBS options
BOUT_walltime         = '06:00:00'
BOUT_run_name         = theRunName
post_process_run_name = 'post' + theRunName.capitalize()
# *****************************************************************************

expandRunner = PBS_runner(\
                directory  = directory ,\
                nproc      = nproc ,\
                # Set spatial domain
                nz         = nz,\
                # Set temporal domain
                nout       = nout  ,\
                timestep   = timestep,\
                # Copy the source file
                cpy_source = True  ,\
                make       = make  ,\
                restart    = restart,\
                restart_from = restart_from,\
                additional = [
                              ('tag',theRunName,0),\
                              ('ownFilters'  , 'type', ownFilterType),\
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

dmp_folders, PBS_ids = expandRunner.execute_runs(\
                         remove_old               = remove_old  ,\
                         post_processing_function = postProcess2,\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
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


# Linear and turb options
# =============================================================================
tSlice              = slice(-200, None, 10)
saveTerms           = False
useHyperViscAzVortD = [True]
timestep            = [1]
BOUT_walltime       = '48:00:00'
# =============================================================================


# The linear runner
# =============================================================================
# Linear options
# *****************************************************************************
includeNoise  = True
forceAddNoise = True
# From previous outputs
restart_from = dmp_folders[0]
# Set the temporal domain
nout = [500]
# Name
theRunName = "a1-KiwiFlat-2-linearPhase1"
# PBS options
BOUT_run_name         = theRunName
post_process_run_name = 'post' + theRunName.capitalize()
# *****************************************************************************

linearRun = PBS_runner(\
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
            additional = [
                          ('tag',theRunName,0),\
                          ('switch'      , 'includeNoise'       , includeNoise ),\
                          ('switch'      , 'forceAddNoise'      ,forceAddNoise),\
                          ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                          ('switch'      , 'saveTerms'          ,saveTerms),\
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

dmp_folders, PBS_ids = linearRun.execute_runs(\
                         remove_old               = remove_old,\
                         post_processing_function = postProcess2,\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
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





# Create the runner
# =============================================================================
# The options for the run
# *****************************************************************************
restart_from = dmp_folders[0]
# Set the temporal domain
nout       = [5000]
# Name
theRunName = "a1-KiwiFlat-3-turbulentPhase1"
# PBS options
BOUT_run_name         = theRunName
post_process_run_name = 'post' + theRunName.capitalize()
# *****************************************************************************

turboRun = PBS_runner(\
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
                additional = [
                              ('tag',theRunName,0),\
                              ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                              ('switch'      , 'saveTerms'          ,saveTerms),\
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

_, _ = turboRun.execute_runs(\
                         remove_old               = remove_old,\
                         post_processing_function = postProcess2,\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
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
