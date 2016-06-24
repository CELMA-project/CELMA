#!/usr/bin/env python

"""Driver which checks the plots."""

from bout_runners import basic_runner
import numpy as np

# The options for the run
# =============================================================================
# Set the temporal domain
restart    = None
remove_old = True
directory  = "data"
# Shall we make?
make       = True
# The number of processors
nproc = 4
nz =1
# =============================================================================


# Create the runner
# =============================================================================
myRuns = basic_runner(\
                      directory  = directory ,\
                      nproc      = nproc     ,\
                      nz = nz,\
                      # Copy the source file
                      cpy_source = True      ,\
                      make       = make      ,\
                      )
# =============================================================================


# Perform the run
# =============================================================================
myRuns.execute_runs(remove_old = remove_old)
# =============================================================================
