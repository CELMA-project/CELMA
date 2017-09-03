#!/usr/bin/env python

"""
Submits refreshDates() to the PBS queue.

**NOTE**: Both refreshDates.py and refreshDatesPBSDriver.py must be
          copied to the root folder for the refresh process.
"""

from common.CELMAPy.driverHelpers import PBSSubmitter
from refreshDates import refreshDates

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
# =============================================================================
sub = PBSSubmitter()
sub.setNodes(nodes=1, ppn=20)
sub.setQueue(queue)
sub.setWalltime("01:00:00")
sub.setMisc(logPath = "." ,\
            mail    = mail,\
            account = account)

sub.setJobName("refreshDates")
sub.submitFunction(refreshDates)
