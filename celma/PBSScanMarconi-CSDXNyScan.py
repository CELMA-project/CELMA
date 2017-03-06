#!/usr/bin/env python

"""
Driver which runs using PBS.
Note: This is done in the old school bout_runners way
"""

from bout_runners import PBS_runner
import numpy as np

mode       = "expand"
saveFolder = mode

# Calculate the ny parameters
ly   = 2.8
rhos = 1.6956e-02 # Copied from logfiles
ratios = np.linspace(0.0975,0.3,5)
# NOTE: ny=66 is the standard
nys  = tuple(int(ny) for ny in ratios*(ly/rhos))

PBSOptions = {\
              "BOUT_walltime" : "72:00:00"    ,\
              "BOUT_queue"    : "xfualongprod",\
              "BOUT_account"  : "FUA11_SOLF"  ,\
              "BOUT_ppn"      : 36            ,\
             }

nproc = 16
if mode == "initial":
    #{{{initial
    # The options for the run
    # ==========================================================================
    # Set the temporal domain
    restart    = None
    remove_old = False
    tempOptions = {\
                   "nout"     : 2  ,\
                   "timestep" : 2e3,\
                  }
    directory  = "CSDXNyScan"
    # Shall we make?
    make       = False
    # ==========================================================================

    for ny in nys:
        while ny%2 != 0:
            ny += 1

        # The PBS options
        # ======================================================================
        # Specify the numbers used for the BOUT runs
        PBSOptions["BOUT_run_name"] = saveFolder + str(ny)
        PBSOptions["BOUT_nodes"] = int(np.ceil(nproc/PBSOptions["BOUT_ppn"]))
        # ======================================================================

        # Create the runner
        # ======================================================================
        myRuns = PBS_runner(\
                    directory  = directory,\
                    nproc      = nproc    ,\
                    # Spatial options
                    ny         = ny       ,\
                    nz         = 1        ,\
                    # Set temporal domain
                    **tempOptions         ,\
                    # Copy the source file
                    cpy_source = True     ,\
                    make       = make     ,\
                    restart    = restart  ,\
                    # PBS options
                    **PBSOptions          ,\
                    # Additional
                    additional = (("tag",saveFolder,0),\
                                  ("ownFilters", "type", "none"),\
                                  ("switch", "useHyperViscAzVortD",False),\
                                 ),\
                    )
        # ======================================================================

        # Perform the run
        # ======================================================================
        myRuns.execute_runs(remove_old = remove_old)
        # ======================================================================
    #}}}
elif mode == "expand":
    #{{{expand
    # The options for the run
    # ==========================================================================
    # Set the temporal domain
    restart    = "overwrite"
    remove_old = False
    tempOptions = {\
                   "nout"     : 2 ,\
                   "timestep" : 25,\
                  }
    directory  = "CSDXNyScan"
    # Shall we make?
    make       = False
    # ==========================================================================

    # The PBS options
    # ==========================================================================
    # Specify the numbers used for the BOUT runs
    PBSOptions["BOUT_run_name"] = saveFolder
    # ==========================================================================

    for ny in nys:
        while ny%2 != 0:
            ny += 1
        restart_from = ("CSDXNyScan/nout_2_timestep_2000.0/ny_{}_nz_1/"
                        "ownFilters_type_none_"
                        "switch_useHyperViscAzVortD_False_tag_initial_0/").\
                       format(ny)

        PBSOptions["BOUT_nodes"] = int(np.ceil(nproc/PBSOptions["BOUT_ppn"]))

        # Create the runner
        # ======================================================================
        myRuns = PBS_runner(\
                    directory  = directory     ,\
                    nproc      = nproc         ,\
                    # Spatial options
                    ny         = ny            ,\
                    nz         = 256           ,\
                    # Set temporal domain
                    **tempOptions              ,\
                    # Copy the source file
                    cpy_source   = True        ,\
                    make         = make        ,\
                    restart      = restart     ,\
                    restart_from = restart_from,\
                    # PBS options
                    **PBSOptions               ,\
                    # Tag (used to catalogize the runs)
                    # Additional
                    additional = (("ownFilters", "type", "none"),\
                                  ("switch", "useHyperViscAzVortD",False),\
                                  ("tag",saveFolder,0),\
                                 ),\
                    )
        # ======================================================================

        # Perform the run
        # ======================================================================
        myRuns.execute_runs(remove_old = remove_old)
        # ======================================================================
    #}}}
