#!/usr/bin/env python

"""
Restarts walltimed run.
NOTE: Copy to correct folder must be done manually.
"""

from bout_runners import PBS_runner
import numpy as np

# The options for the run
# =============================================================================
# *****************************************************************************
B0 = (  8.0e-2,)
Lx = (  6.2906,)
Ly = (220.1715,)
# *****************************************************************************
nz = 1
restart      = "append"
restart_from = "CSDXMagFieldScanAr/nout_2_timestep_2000.0/nz_1/geom_Lx_6.2906_geom_Ly_220.1715_input_B0_0.08_ownFilters_type_none_switch_useHyperViscAzVortD_False_tag_CSDXMagFieldScanAr-0-initialize_0/"
timestep     = (2e3,)
nout         = (1,)
directory    = "CSDXMagFieldScanAr"
theRunName   = "restart-{}-0-initialize".format(directory)
make         = False
# *****************************************************************************
ownFilterType       = "none"
useHyperViscAzVortD = False
# *****************************************************************************
# =============================================================================

# The PBS options
# =============================================================================
# Specify the numbers used for the BOUT runs
nproc                 = 48
BOUT_nodes            = 3
BOUT_ppn              = 20
BOUT_walltime         = "48:00:00"
BOUT_run_name         = theRunName
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
            cpy_source = True  ,\
            make       = make  ,\
            restart    = restart,\
            restart_from = restart_from,\
            series_add   = (\
                ("input", "B0", B0),\
                ("geom" , "Lx", Lx),\
                ("geom" , "Ly", Ly),\
                           ),\
            additional   = (\
                ("tag"       , theRunName           , 0),\
                ("ownFilters", "type"               , ownFilterType),\
                ("switch"    , "useHyperViscAzVortD", useHyperViscAzVortD),\
                           ),\
            # PBS options
            BOUT_nodes    = BOUT_nodes           ,\
            BOUT_ppn      = BOUT_ppn             ,\
            BOUT_walltime = BOUT_walltime        ,\
            BOUT_run_name = BOUT_run_name        ,\
            )
# =============================================================================


# Perform the run
# =============================================================================
myRuns.execute_runs()
# =============================================================================
