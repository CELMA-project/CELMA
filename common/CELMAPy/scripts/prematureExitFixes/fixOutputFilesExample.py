#!/usr/bin/env python

"""
Example on how to fix output files of premature exit runs.
Run this in case some of the processors wrote a timestep, others didn't.

**NOTE**: DO NOT RUN THIS unless tere is a problem.
"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.driverHelpers import PBSSubmitter
from CELMAPy.repairBrokenExit import repairBrokenExit

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
# This is the path for the files where the output is broken
path = ("CSDXMagFieldScanAr/nout_5000_timestep_1/"
        "geom_Lx_3.1453_geom_Ly_110.0858_input_B0_0.04_"
        "switch_saveTerms_False_switch_useHyperViscAzVortD_"
        "True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_7/")
# =============================================================================

directory = path.split("/")[0]

sub = PBSSubmitter()
sub.setNodes(nodes=1, ppn=20)
sub.setQueue(queue)
sub.setWalltime("01:00:00")
sub.setMisc(logPath = os.path.join(directory,"postLogs"),\
            mail    = mail,\
            account = account)

sub.setJobName("repairBrokenExit")
sub.submitFunction(repairBrokenExit, args=(path,), kwargs={})
