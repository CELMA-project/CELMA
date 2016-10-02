#!/usr/bin/env python

"""Method of exact solution driver"""

from bout_runners.bout_runners import basic_runner

# The options for the run
# =============================================================================
# Additional options
remove_old = True
make       = False
nproc      = 4
nn         = [1.0e15, 5.0e15, 1.0e16, 5.0e16, 1.0e17]
# =============================================================================


# Create the runner
# =============================================================================
my_runs = basic_runner(\
            nproc      = nproc ,\
            # Copy the source file
            additional = ('input', 'nn', nn),\
            cpy_source = True  ,\
            make       = make  ,\
            )
# =============================================================================


# Perform the run
# =============================================================================
my_runs.execute_runs(remove_old = remove_old)
# =============================================================================
