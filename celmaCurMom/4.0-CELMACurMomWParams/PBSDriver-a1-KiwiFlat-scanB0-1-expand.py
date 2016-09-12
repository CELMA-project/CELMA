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
ownFilterType = "none"
B0s = [9.0e-2  , 8.0e-2   , 7.0e-2 , 6.0e-2  , 5.0e-2  ]
Lxs = [4.3466  , 3.8637   , 3.3807 , 2.8978  , 2.4148  ]
Lys = [243.4121, 216.3663, 189.3205, 162.2747, 135.2289]
# *****************************************************************************
# Set the spatial domain
nz = 256
# Set the temporal domain
remove_old = False
restart    = "overwrite"
# Uncomment this if you just want to plot
# restart      = None;
# Set the temporal domain
timestep   = [5]
nout       = [2]
directory  = "a1-KiwiFlat"
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
tSlice     = slice(-2, None)
showPlot   = False
savePlot   = True
# =============================================================================


# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 48
BOUT_nodes            = 3
BOUT_ppn              = 16
BOUT_walltime         = '06:00:00'
post_process_nproc    = 1
post_process_nodes    = 1
post_process_ppn      = 20
post_process_walltime = '0:29:00'
post_process_queue    = 'xpresq'
# =============================================================================


for Lx, Ly, B0 in zip(Lxs, Lys, B0s):
    restart_from  = "a1-KiwiFlat/nout_2_timestep_2000.0/nz_1/geom_Lx_{}_geom_Ly_{}_input_B0_{}_ownFilters_type_none_tag_a1-KiwiFlat-scanB0-0-initialize_0/".format(Lx, Ly, B0)
    theRunName    = "a1-KiwiFlat-B0-{}-1-expand".format(B0)
    BOUT_run_name = theRunName
    post_process_run_name = 'post' + theRunName.capitalize()
    # Create the runner
    # =========================================================================
    myRuns = PBS_runner(\
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
                series_add = [
                              ('input', 'B0', B0),\
                              ('geom', 'Lx', Lx),\
                              ('geom', 'Ly', Ly),\
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
    # =========================================================================


    # Perform the run
    # =========================================================================
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
                         tSlice         = tSlice            ,\
                         savePlot       = savePlot          ,\
                         saveFolderFunc = "scanWTagSaveFunc",\
                         theRunName     = theRunName        ,\
                        )
    # =========================================================================
