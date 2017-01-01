#!/usr/bin/env python

"""Driver which runs using PBS."""

import pickle
from time import sleep

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../common")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.driverHelpers import PBSSubmitter, pathMerger
from standardPlots import fourierModesPlot

#{{{runFourierModes
def runFourierModes():
    """
    Runs the fourier modes
    """
    loopOver = zip(dmpFolders["turbulence"], paramKeys, rJobs)
    for dmp_folders, key, nr in loopOver:

        collectPaths = mergeFromLinear[key]
        dmp_folders   = (dmp_folders,)
        # fourierModesPlot(dmp_folders, collectPaths)

        args = (dmp_folders, collectPaths)
        sub.setJobName("fourierModes{}".format(nr))
        sub.submitFunction(fourierModesPlot, args=args)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{runFourierModesSliced
def runFourierModesSliced():
    loopOver = zip(dmpFolders["turbulence"], paramKeys, rJobs)
    for dmp_folders, key, nr in loopOver:

        # Find tSlice
        found = False
        for tkey in tSlices.keys():
            if tkey in dmp_folders:
                tSlice = tSlices[tkey]
                found = True
                break
        if not(found):
            raise ValueError("Could not find correct slice")

        collectPaths = mergeFromLinear[key]
        dmp_folders   = (dmp_folders,)
        args = (dmp_folders, collectPaths)
        kwargs = {"tSlice":(tSlice,)}
        sub.setJobName("fourierModesSliced{}".format(nr))
        sub.submitFunction(fourierModesPlot, args=args, kwargs=kwargs)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{Globals
directory = "CSDXMagFieldScanAr"

with open(os.path.join(directory, "dmpFoldersDict.pickle"), "rb") as f:
        dmpFolders = pickle.load(f)

mergeFromLinear =\
        pathMerger(dmpFolders, ("linear", "turbulence", "extraTurbulence"))

paramKeys = tuple(sorted(list(mergeFromLinear.keys())))
rJobs     = range(len(paramKeys))

# Generate the submitter
sub = PBSSubmitter()
sub.setNodes(nodes=1, ppn=2)
sub.setQueue("xpresq")
sub.setWalltime("00:05:00")
#}}}

if __name__ == "__main__":

    # Run the fourier modes
    runFourierModes()
    # Set linear slices and plot the sliced fourier modes
    tSlices = {\
               "B0_0.02":slice(80,240)  ,\
               "B0_0.04":slice(800,1250),\
               "B0_0.06":slice(180,300) ,\
               "B0_0.08":slice(100,225) ,\
               "B0_0.1" :slice(80,210)  ,\
               }
    runFourierModesSliced()
