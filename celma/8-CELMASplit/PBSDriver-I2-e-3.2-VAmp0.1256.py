#!/usr/bin/env python

"""Driver which runs using PBS."""

from bout_runners import PBS_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common/python')
# Sys path is a list of system paths
sys.path.append(commonDir)

from plotting import combinedDriver

# The options for the run
# =============================================================================
# *****************************************************************************
ownOpType           = "BasicBrackets"
ownFilterType       = "none"
sAmp                = [0.020]
saveDdt             = True
includeNoise        = False
forceAddNoise       = False
useHyperViscAzVortD = True
artHyperAzVortD     = [1e-1]
# *****************************************************************************
remove_old = False
restart    = "overwrite"
# Uncomment this if you just want to plot
# restart      = None;
restart_from = "f-moreSource3Brackets/nout_300_timestep_10/nz_256/cst_artHyperAzVortD_0.1_ownFilters_type_none_ownOperators_type_BasicBrackets_switch_forceAddNoise_True_switch_includeNoise_True_switch_saveDdt_True_switch_useHyperViscAzVortD_True_tag_2-e-3.2-moreSourceAddNoiseHyperviscBasicBrackets_0_theSource_a_0.02/"
# Set the spatial domain
nz = 256
# Set the temporal domain
nout       = [100]
timestep   = [1]
directory  = "f-moreSource3Brackets"
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
theRunName = "I2-e-3.2-VAmp0.1256"
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
                          ('ownOperators', 'type', ownOpType    ),\
                          ('ownFilters'  , 'type', ownFilterType),\
                          ('theSource', 'a', sAmp               ),\
                          ('switch'      , 'saveDdt'            , saveDdt      ),\
                          ('switch'      , 'includeNoise'       , includeNoise ),\
                          ('switch'      , 'forceAddNoise'      ,forceAddNoise),\
                          ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
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
                     savePlot       = savePlot          ,\
                     saveFolderFunc = "scanWTagSaveFunc",\
                     theRunName     = theRunName        ,\
                    )
# =============================================================================
