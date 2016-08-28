#!/usr/bin/env python

"""Method of exact solution driver"""

from bout_runners.bout_runners import basic_runner

# The options for the run
# =============================================================================
# Additional options
remove_old = True
make       = False
nproc      = 4
# =============================================================================


# Create the runner
# =============================================================================
my_runs = basic_runner(\
            nproc      = nproc ,\
            # Copy the source file
            cpy_source = True  ,\
            make       = make  ,\
            )
# =============================================================================


# Perform the run
# =============================================================================
my_runs.execute_runs(remove_old = remove_old)
# =============================================================================
