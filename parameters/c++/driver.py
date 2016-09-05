#!/usr/bin/env python

"""Method of exact solution driver"""

from bout_runners.bout_runners import basic_runner

# The options for the run
# =============================================================================
# Additional options
remove_old = True
make       = False
nproc      = 4
B0 = [9e-2, 8e-2, 7e-2, 6e-2, 5e-2]
# =============================================================================


# Create the runner
# =============================================================================
my_runs = basic_runner(\
            nproc      = nproc ,\
            # Copy the source file
            additional = ('input', 'B0', B0),\
            cpy_source = True  ,\
            make       = make  ,\
            )
# =============================================================================


# Perform the run
# =============================================================================
my_runs.execute_runs(remove_old = remove_old)
# =============================================================================
