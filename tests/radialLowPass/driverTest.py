#!/usr/bin/env python

"""Test driver"""

from bout_runners import basic_runner
import numpy as np

from postProcess.checkModes import checkModes as postProcess

# The options for the run
# =============================================================================
# Additional options
remove_old = True
make       = True
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
my_runs.execute_runs(\
                     remove_old = remove_old,\
                     # Set the proper directory
                     post_processing_function = postProcess,\
                    )
# =============================================================================
