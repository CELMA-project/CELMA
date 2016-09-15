#!/usr/bin/env python

"""Driver which checks the plots."""

from bout_runners import basic_runner
import numpy as np

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath('./../common')
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPython.drivers import postBoutRunner

# If you just want to post-process
justPostProcess = True

# The options for the run
# =============================================================================
# Set the temporal domain
restart    = None
remove_old = False
nout       = 2
timestep   = 1e-10
directory  = "a1-KiwiFlat"
# Shall we make?
make       = False
# The number of processors
nproc = 4
# =============================================================================

collectionFolders = []

# Create the runner and run
# =============================================================================
myRuns = basic_runner(\
                      directory  = directory ,\
                      nproc      = nproc     ,\
                      # Set temporal domain
                      nout       = nout      ,\
                      timestep   = timestep  ,\
                      # Copy the source file
                      cpy_source = True      ,\
                      make       = make      ,\
                      restart    = restart   ,\
                      )

dmp_folder, _ = myRuns.execute_runs(remove_old = remove_old)
collectionFolders.append(dmp_folder[0])
# =============================================================================
if justPostProcess:
    restart = None
else:
    restart    = "overwrite"

# The options for the run
# =============================================================================
# Set the temporal domain
nout       = 3
# =============================================================================

# Create the runner and run
# =============================================================================
myRuns = basic_runner(\
                      directory  = directory ,\
                      nproc      = nproc     ,\
                      # Set temporal domain
                      nout       = nout      ,\
                      timestep   = timestep  ,\
                      # Copy the source file
                      cpy_source = True      ,\
                      make       = make      ,\
                      restart    = restart   ,\
                      restart_from = dmp_folder[0]\
                      )

dmp_folder, _ = myRuns.execute_runs(remove_old = remove_old)
collectionFolders.append(dmp_folder[0])
# =============================================================================

# The options for the run
# =============================================================================
# Set the temporal domain
nout       = 10
# =============================================================================


# Create the runner and run
# =============================================================================
myRuns = basic_runner(\
                      directory  = directory ,\
                      nproc      = nproc     ,\
                      # Set temporal domain
                      nout       = nout      ,\
                      timestep   = timestep  ,\
                      # Copy the source file
                      cpy_source = True      ,\
                      make       = make      ,\
                      restart    = restart   ,\
                      restart_from = dmp_folder[0],\
                      # Additional
                      additional = [\
                        ('switch', 'includeNoise' , True),\
                        ('switch', 'forceAddNoise', True),\
                                   ]
                      )

dmp_folder, _ = myRuns.execute_runs(remove_old = remove_old)
collectionFolders.append(dmp_folder[0])
# =============================================================================


# Plot the probe data
postBoutRunner(dmp_folder,\
               # postBoutRunner input
               driverName = "plotProbes",\
               # PostProcessDriver input
               convertToPhysical = False             ,\
               # subPolAvg         = False             ,\
               showPlot          = False             ,\
               savePlot          = True              ,\
               # saveFolder        = None              ,\
               saveFolderFunc    = "scanWTagSaveFunc",\
               # useSubProcess     = True              ,\
               theRunName        = "probeTest"       ,\
               # StatsAndSignalsDrivers input
               paths             = collectionFolders,\
               # DriversProbes input
               var               = 'n'               ,\
               yInd              = 2*8               ,\
               nProbes           = 5                 ,\
               # Choose this in order to have some gradients
               steadyStatePath   = collectionFolders[1] ,\
               maxMode           = 7                 ,\
              )
