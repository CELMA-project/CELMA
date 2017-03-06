#!/usr/bin/env python

"""
Example on how to fix broken exit runs.
In this case some of the processors wrote a timestep, others didn't.

**NOTE**: DO NOT RUN THIS unless tere is a problem.
"""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.driverHelpers import PBSSubmitter
from CELMAPy.repairBrokenExit import repairBrokenExit

directory = "CSDXNeutralScanAr"

sub = PBSSubmitter()
sub.setNodes(nodes=1, ppn=20)
sub.setQueue("xfualongprod")
sub.setWalltime("01:00:00")
sub.setMisc(logPath = os.path.join(directory,"postLogs"),\
            mail    = "mmag@fysik.dtu.dk",\
            account = "FUA11_SOLF")

path = ("CSDXMagFieldScanAr/nout_5000_timestep_1/"
        "geom_Lx_3.1453_geom_Ly_110.0858_input_B0_0.04_"
        "switch_saveTerms_False_switch_useHyperViscAzVortD_"
        "True_tag_CSDXMagFieldScanAr-3-turbulentPhase1_0/restart_7/")

sub.setJobName("repairBrokenExit")
sub.submitFunction(repairBrokenExit, args=(path,), kwargs={})
