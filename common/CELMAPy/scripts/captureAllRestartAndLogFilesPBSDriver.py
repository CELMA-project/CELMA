#!/usr/bin/env python

"""
Submits captureAllRestartAndLogFiles() to the PBS queue.

**NOTE**: Both captureAllRestartAndLogFiles.py and
          captureAllRestartAndLogFilesPBSDriver.py must be copied to the
          desired folder.
"""

from common.CELMAPy.driverHelpers import PBSSubmitter
from captureAllRestartAndLogFiles import captureAllRestartAndLogFiles

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
# Directory to recursively find log and restart files from.
directory = "."
# =============================================================================
sub = PBSSubmitter()
sub.setNodes(nodes=1, ppn=20)
sub.setQueue(queue)
sub.setWalltime("01:00:00")
sub.setMisc(logPath = "." ,\
            mail    = mail,\
            account = account)

sub.setJobName("captureAllRestartAndLogFiles")
sub.submitFunction(captureAllRestartAndLogFiles, args = (directory,))
