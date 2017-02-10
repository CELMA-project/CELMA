#!/usr/bin/env python

"""
Driver which runs using PBS.
Note: This is done in the old school bout_runners way
"""

from bout_runners import PBS_runner

mode       = "initial"
saveFolder = mode


# Calculate the ny parameters
ly   = 2.8
rhos = 1.6956e-02 # Copied from logfiles
# NOTE: ny=66 is the standard
nys  = tuple(int(ny) for ny in a*(ly/rhos))

if mode == "initial":
    #{{{initial
    # The options for the run
    # ==========================================================================
    # Set the temporal domain
    restart    = None
    remove_old = False
    nout       = [2]
    timestep   = [2e3]
    directory  = "CSDXNyScan"
    # Shall we make?
    make       = False
    # ==========================================================================

    # The PBS options
    # ==========================================================================
    # Specify the numbers used for the BOUT runs
# FIXME:
    nproc                 = 96
    BOUT_nodes            = 1
    BOUT_ppn              = 20
    BOUT_walltime         = "72:00:00"
    BOUT_queue            = "xfualong"
    BOUT_run_name         = saveFolder
    # ==========================================================================

    # Create the runner
    # ==========================================================================
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
                # PBS options
                BOUT_nodes            = BOUT_nodes   ,\
                BOUT_ppn              = BOUT_ppn     ,\
                BOUT_queue            = BOUT_queue   ,\
                BOUT_walltime         = BOUT_walltime,\
                BOUT_run_name         = BOUT_run_name,\
                # Tag (used to catalogize the runs)
                additional = ("tag",saveFolder,0)    ,\
                )
    # ==========================================================================

    # Perform the run
    # ==========================================================================
    myRuns.execute_runs(remove_old = remove_old)
    # ==========================================================================
    #}}}
elif mode == "steadyState":
    #{{{steadyState
    raise NotImplementedError("steadyState is to be implemented")
    #}}}
