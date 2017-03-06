#!/usr/bin/env python

"""Driver which plots the results of the simulations."""

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.driverHelpers import PBSSubmitter
from standardPlots import fields1DAnimation

# INPUT
# =============================================================================
# If you would like a mail on finished job enter your mail here
# Example: "john@doe.com"
mail  = None
# If the queuing system uses accounts, set the account here
# Example: "FUA11_SOLF"
account = None
# Usually, the queueing system has its own default queue, if not,
# specify here
# Example: "xfualongprod"
queue = None
# =============================================================================

directory = "CSDXNyScan"

# Generate the submitter
sub = PBSSubmitter()
sub.setNodes(nodes=1, ppn=20)
sub.setQueue(queue)
sub.setWalltime("00:15:00")
sub.setMisc(logPath = os.path.join(directory,"postLogs"),\
            mail    = mail,\
            account = account)
sub.setQueue(queue)
sub.toggleSubmitOrRun()

# Define plotting parameters
plotSuperKwargs = {\
                   "extension"       : "pdf"     ,\
                   "showPlot"        : False     ,\
                   "savePlot"        : True      ,\
                   "savePath"        : None      ,\
                   "savePathFunc"    : "onlyScan",\
                   "timeStampFolder" : False     ,\
                   "scanParameter"   : "ny"      ,\
                  }

nys= (16, 24, 32, 42, 50)
dmpFolders = tuple((directory+"/nout_2_timestep_25/ny_{}_nz_256/"
              "ownFilters_type_none_"
              "switch_useHyperViscAzVortD_False_tag_expand_0/").\
              format(ny) for ny in nys)

for nr, dmp_folders in enumerate(dmpFolders):
    dmp_folders  = (dmp_folders,)
    args = (dmp_folders, dmp_folders, plotSuperKwargs)
    kwargs = {"hyperIncluded" : False,\
              "boussinesq"    : False,\
              "tSlice"        : slice(-1,-1),\
              "yInd"          : 4,\
             }
    sub.setJobName("steadyState{}".format(nr))
    sub.submitFunction(fields1DAnimation, args=args, kwargs=kwargs)
