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

# If you just want to post-process
justPostProcess = False

#{{{restart_from_func
def restart_from_func(dmp_folder, **kwargs):
    """
    Function which returns the restart from folder

    Parameters
    ----------
    dmp_folder : str
        Given by the bout_runners
    one_of_the_restart_paths_in_scan : str
        One of the restart paths from a previously run scan. This
        paramemter will be given as a kwargs
    kwargs : dict
        Dictionary with additional keyword arguments, given by
        bout_runners.
        One of the arguments (given as kwargs to execute_runs) is
        one_of_the_restart_paths_in_scan.
    """

    # Ensure that the variable is set
    one_of_the_restart_paths_in_scan =\
            kwargs["one_of_the_restart_paths_in_scan"]

    # Values in the current scan
    valNames = ["Lx", "Ly", "B0"]

    # Make a template of the restart paths
    # Split the dmp_folder at _ (separator between the variable and its
    # value)
    restart_template = one_of_the_restart_paths_in_scan.split("_")

    # Remove the values from one_of_the_restart_paths_in_scan in order
    # to make it possible to enter values using the str.format()
    # function
    for valName in valNames:
        indices = [ind for ind, el in enumerate(restart_template)\
                   if el == valName]

        for ind in indices:
            # The value is found to the right of the variable
            restart_template[ind + 1] = "{{0[{}]}}".format(valName)

    # Join the restart_template to one string
    restart_template = "_".join(restart_template)

    # Find the values to put into the template
    splitted = dmp_folder.split("_")
    values = {}
    for valName in valNames:
        # Find the first index of val
        ind = splitted.index(valName)
        # The value is found to the right of the variable
        values[valName] = splitted[ind+1]

    # Insert the values into the restart_template string using the
    # values dictionary
    restart_from = restart_template.format(values)

    return restart_from
#}}}

# Common options
# =============================================================================
# The scan
# *****************************************************************************
B0 = [1.0e-1  , 9.0e-2  , 8.0e-2   , 7.0e-2 , 6.0e-2  , 5.0e-2  ]
Lx = [4.8296  , 4.3466  , 3.8637   , 3.3807 , 2.8978  , 2.4148  ]
Ly = [270.4579, 243.4121, 216.3663, 189.3205, 162.2747, 135.2289]
# *****************************************************************************

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
directory  = "a1-KiwiFlat"
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

dmp_folders, PBS_ids = initRunner.execute_runs(\
                         remove_old               = remove_old,\
                         post_processing_function = postBoutRunner,\
                         # This function will be called every time after
                         # performing a run
                         post_process_after_every_run = True,\
                         # Below are the kwargs arguments being passed to
                         # the post processing function
                         # Switches
                         driverName     = "plot1D2DAndFluctDriver",\
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

# Expand, Linear and turb options
# =============================================================================
if justPostProcess:
    restart = None
else:
    restart = "overwrite"
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
one_of_the_restart_paths_in_scan = dmp_folders[0]
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
                restart_from = restart_from_func,\
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

dmp_folders, PBS_ids = expandRunner.execute_runs(\
                         remove_old               = remove_old  ,\
                         post_processing_function = postBoutRunner,\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
                         # This function will be called every time after
                         # performing a run
                         post_process_after_every_run = True,\
                         # Below are the kwargs arguments being passed to
                         # the post processing function
                         # Switches
                         driverName     = "plot1D2DAndFluctDriver",\
                         xguards        = xguards           ,\
                         yguards        = yguards           ,\
                         xSlice         = xSlice            ,\
                         ySlice         = ySlice            ,\
                         zSlice         = zSlice            ,\
                         tSlice         = tSlice            ,\
                         savePlot       = savePlot          ,\
                         saveFolderFunc = "scanWTagSaveFunc",\
                         theRunName     = theRunName        ,\
                         # Below are the kwargs given to the
                         # restart_from_func
                         one_of_the_restart_paths_in_scan =\
                         one_of_the_restart_paths_in_scan,\
                        )
# =============================================================================


# Linear and turb options
# =============================================================================
tSlice              = slice(-200, None, 10)
saveTerms           = False
useHyperViscAzVortD = [True]
timestep            = [1]
BOUT_walltime       = '72:00:00'
# =============================================================================


# The linear runner
# =============================================================================
# Linear options
# *****************************************************************************
includeNoise  = True
forceAddNoise = True
# From previous outputs
one_of_the_restart_paths_in_scan = dmp_folders[0]
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
            restart_from = restart_from_func,\
            additional = [
                          ('tag',theRunName,0),\
                          ('switch'      , 'includeNoise'       , includeNoise ),\
                          ('switch'      , 'forceAddNoise'      ,forceAddNoise),\
                          ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                          ('switch'      , 'saveTerms'          ,saveTerms),\
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

dmp_folders, PBS_ids = linearRun.execute_runs(\
                         remove_old               = remove_old,\
                         post_processing_function = postBoutRunner,\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
                         # This function will be called every time after
                         # performing a run
                         post_process_after_every_run = True,\
                         # Below are the kwargs arguments being passed to
                         # the post processing function
                         # Switches
                         driverName     = "plot1D2DAndFluctDriver",\
                         xguards        = xguards           ,\
                         yguards        = yguards           ,\
                         xSlice         = xSlice            ,\
                         ySlice         = ySlice            ,\
                         zSlice         = zSlice            ,\
                         tSlice         = tSlice            ,\
                         savePlot       = savePlot          ,\
                         saveFolderFunc = "scanWTagSaveFunc",\
                         theRunName     = theRunName        ,\
                         # Below are the kwargs given to the
                         # restart_from_func
                         one_of_the_restart_paths_in_scan =\
                         one_of_the_restart_paths_in_scan,\
                        )
# =============================================================================





# Create the runner
# =============================================================================
# The options for the run
# *****************************************************************************
one_of_the_restart_paths_in_scan = dmp_folders[0]
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
                restart_from = restart_from_func,\
                additional = [
                              ('tag',theRunName,0),\
                              ('switch'      , 'useHyperViscAzVortD',useHyperViscAzVortD),\
                              ('switch'      , 'saveTerms'          ,saveTerms),\
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

_, _ = turboRun.execute_runs(\
                         remove_old               = remove_old,\
                         post_processing_function = postBoutRunner,\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
                         # This function will be called every time after
                         # performing a run
                         post_process_after_every_run = True,\
                         # Below are the kwargs arguments being passed to
                         # the post processing function
                         # Switches
                         driverName     = "plot1D2DAndFluctDriver",\
                         xguards        = xguards           ,\
                         yguards        = yguards           ,\
                         xSlice         = xSlice            ,\
                         ySlice         = ySlice            ,\
                         zSlice         = zSlice            ,\
                         tSlice         = tSlice            ,\
                         savePlot       = savePlot          ,\
                         saveFolderFunc = "scanWTagSaveFunc",\
                         theRunName     = theRunName        ,\
                         # Below are the kwargs given to the
                         # restart_from_func
                         one_of_the_restart_paths_in_scan =\
                         one_of_the_restart_paths_in_scan,\
                        )
# =============================================================================
