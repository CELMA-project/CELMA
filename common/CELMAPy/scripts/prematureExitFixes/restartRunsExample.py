#!/usr/bin/env python

"""
Example driver which continues run after premature exit.

In this example the cluster terminated prematurely in the initialization
phase.

This script uses the native bout_runners.

To fix the logfiles and restart files, run
fixLogFilesExample.py and
fixRestartFilesExample.py
when the runs started here has finished.

**NOTE**: The code must be altered if the premature exit occured during
          the expand phase.
          If the crash happened after the expand phase, the keyword
          argument
          "restartTurb = desiredRestarts"
          of "ScanDriverObject.runScan()" may be used instead.
"""

from bout_runners import PBS_runner

# INPUT
# =============================================================================
# If the queuing system uses accounts, set the account here
# Example: "FUA11_SOLF"
account = None
# Usually, the queueing system has its own default queue, if not,
# specify here
# Example: "xfualongprod"
queue = None
# If you would like a mail on finished job enter your mail here
# Example: "john@doe.com"
mail = None
# Restart options
directory = "CSDXMagFieldScanAr"
# How long did the SimTime get
reached  = 2000
# What was the desired SimTime
end      = 4000
# How many outputs would you like in the restart runs
nout     = 1000
# Set the scan
B0s = (  1.0e-1,   8.0e-2,   6.0e-2)
Lxs = (  7.8633,   6.2906,   4.7180)
Lys = (275.2144, 220.1715, 165.1286)
# =============================================================================

for B0, Lx, Ly in zip(B0s, Lxs, Lys):
    series_add = (\
                  ("input", "B0", B0),\
                  ("geom" , "Lx", Lx),\
                  ("geom" , "Ly", Ly),\
                 )

    # Set common runner options
    commonRunnerOptions = {\
                           "directory"    : directory,\
                           "nproc"        : 48     ,\
                           "cpy_source"   : True   ,\
                           "BOUT_nodes"   : 2      ,\
                           "BOUT_ppn"     : 32     ,\
                           "BOUT_queue"   : queue  ,\
                           "BOUT_account" : account,\
                           "BOUT_mail"    : mail   ,\
                          }

    timestep = (end-reached)/nout

    #{{{Init options
    # Name
    theRunName = directory + "-0-initialize"
    # Set the spatial domain
    nz = 1
    # Set the temporal domain
    restart = "overwrite"
    restart_from = "CSDXMagFieldScanAr/nout_2_timestep_2000.0/nz_1/geom_Lx_{}_geom_Ly_{}_input_B0_{}_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0BAK".format(Lx[0],Ly[0],B0[0])
    # Filter
    ownFilterType = "none"
    hyper = ("switch", "useHyperViscAzVortD",False)
    # Init options
    initOptions =\
      {\
       "timestep" : (timestep,),\
       "nout"     : (nout,),\
      }
    # PBS options
    initPBSOptions = {\
                        "BOUT_run_name": theRunName ,\
                        "BOUT_walltime" : "72:00:00",\
                     }
    #}}}
    #{{{Run
    initRunner = PBS_runner(\
        make       = False   ,\
        nz         = nz      ,\
        restart    = restart ,\
        restart_from = restart_from,\
        **initOptions        ,\
        # Set additional option
        additional = (
            ("tag",        theRunName, 0            ),\
            ("ownFilters", "type"    , ownFilterType),\
            hyper                                    ,\
                     )               ,\
        series_add = series_add,\
        # PBS options
        **initPBSOptions,\
        # Common options
        **commonRunnerOptions,\
                    )

    init_dmp_folders, init_PBS_ids = initRunner.execute_runs()
    #}}}
