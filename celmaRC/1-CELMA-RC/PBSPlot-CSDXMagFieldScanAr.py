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
from standardPlots import (fourierModesPlot,\
                           growthRatesPlot ,\
                           energyPlot      ,\
                           posOfFluctPlot  ,\
                           PSD2DPlot       ,\
                           skewKurtPlot    ,\
                           )

#{{{findSlices
def findSlices(dmp_folders):
    """
    Finds the slices

    Parameters
    ----------
    dmp_folders : str
        Dump folder under investigation
    """
    found = False
    for tkey in tSlices.keys():
        if tkey in dmp_folders:
            tSlice = tSlices[tkey]
            found = True
            break
    if not(found):
        raise ValueError("Could not find correct slice")

    return tSlice
#}}}

#{{{runFourierModes
def runFourierModes(sliced = False):
    """
    Runs the fourier modes

    Parameters
    ----------
    sliced : bool
        Whether or not to slice the time
    """

    loopOver = zip(dmpFolders["turbulence"],\
                   dmpFolders["expand"],\
                   paramKeys,\
                   rJobs)
    for dmp_folders, steadyStatePath, key, nr in loopOver:

        if sliced:
            tSlice = findSlices(dmp_folders)
            if tSlice is None:
                continue

        collectPaths = mergeFromLinear[key]
        dmp_folders   = (dmp_folders,)
        args = (dmp_folders, collectPaths, steadyStatePath)
        if sliced:
            kwargs = {"tSlice":(tSlice,)}
            sub.setJobName("fourierModesSliced{}".format(nr))
        else:
            kwargs = {}
            sub.setJobName("fourierModes{}".format(nr))
        sub.submitFunction(fourierModesPlot, args=args, kwargs=kwargs)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{runGrowthRates
def runGrowthRates():
    """
    Runs the growth rates
    """

    # NOTE: The ordering of param is in descending order (because of the
    #       organization in PBSScan)
    dmp_folders = (mergeFromLinear["param0"][-1],)
    keys = tuple(sorted(list(mergeFromLinear.keys())))
    scanCollectPaths = tuple(mergeFromLinear[key] for key in keys)

    # NOTE: The ordering is of the keys are in descending order (because
    #       of the organization in PBSScan)
    steadyStatePaths = dmpFolders["expand"]

    # Set the indices
    tKeys     = tuple(sorted(list(tSlices.keys()), reverse=True))
    startInds = tuple(tSlices[key].start for key in tKeys)
    endInds   = tuple(tSlices[key].stop for key in tKeys)

    args = (dmp_folders, scanCollectPaths, steadyStatePaths, scanParameter, startInds, endInds)
    sub.setJobName("growthRates")
    sub.submitFunction(growthRatesPlot, args=args)
#}}}

#{{{runEnergy
def runEnergy(sliced=False):
    """
    Runs the energies

    Parameters
    ----------
    sliced : bool
        Whether or not to slice the time
    """
    loopOver = zip(dmpFolders["turbulence"], paramKeys, rJobs)
    for dmp_folders, key, nr in loopOver:

        if sliced:
            # Find tSlice
            tSlice = findSlices(dmp_folders)
            if tSlice is None:
                continue

        collectPaths = mergeFromLinear[key]
        dmp_folders  = (dmp_folders,)
        args = (dmp_folders, collectPaths)
        if sliced:
            kwargs = {"tSlice":tSlice}
            sub.setJobName("energySliced{}".format(nr))
        else:
            kwargs = {}
            sub.setJobName("energy{}".format(nr))
        sub.submitFunction(energyPlot, args=args, kwargs=kwargs)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{runPosOfFluct
def runPosOfFluct():
    """
    Runs the position of fluct
    """

    loopOver = zip(dmpFolders["turbulence"],\
                   dmpFolders["expand"],\
                   paramKeys,\
                   rJobs)
    for dmp_folders, steadyStatePath, key, nr in loopOver:

        # Find tSlice
        tSlice = findSlices(dmp_folders)
        if tSlice is None:
            continue

        collectPaths = mergeFromLinear[key]
        dmp_folders  = (dmp_folders,)
        args = (dmp_folders, collectPaths, steadyStatePath)
        kwargs = {"tSlice":tSlice}
        sub.setJobName("posOfFluctSliced{}".format(nr))
        sub.submitFunction(posOfFluctPlot, args=args, kwargs=kwargs)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{runPSD2D
def runPSD2D():
    """
    Runs the PSD2D
    """

    loopOver = zip(dmpFolders["turbulence"],\
                   paramKeys,\
                   rJobs)
    for dmp_folders, key, nr in loopOver:

        # Find tSlice
        tSlice = findSlices(dmp_folders)
        if tSlice is None:
            continue

        collectPaths = mergeFromLinear[key]
        dmp_folders  = (dmp_folders,)
        args = (dmp_folders, collectPaths)
        kwargs = {"tSlice":tSlice}
        sub.setJobName("PSD2DPlotSliced{}".format(nr))
        sub.submitFunction(PSD2DPlot, args=args, kwargs=kwargs)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{runSkewKurt
def runSkewKurt():
    """
    Runs the skewness and kurtosis
    """

    loopOver = zip(dmpFolders["turbulence"],\
                   paramKeys,\
                   rJobs)
    for dmp_folders, key, nr in loopOver:

        # Find tSlice
        tSlice = findSlices(dmp_folders)
        if tSlice is None:
            continue

        collectPaths = mergeFromLinear[key]
        dmp_folders  = (dmp_folders,)
        args = (dmp_folders, collectPaths)
        kwargs = {"tSlice":tSlice}
        sub.setJobName("skewnessKurtosisSliced{}".format(nr))
        sub.submitFunction(skewKurtPlot, args=args, kwargs=kwargs)
        # Sleep 1.5 seconds to ensure that tmp files will have different names
        sleep(1.5)
#}}}

#{{{Globals
directory = "CSDXMagFieldScanAr"
scanParameter = "B0"

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
    # runFourierModes(sliced=False)
    # Set linear slices and plot the sliced fourier modes
    tSlices = {\
               "B0_0.02":slice(80,240)  ,\
               "B0_0.04":slice(800,1250),\
               "B0_0.06":slice(180,300) ,\
               "B0_0.08":slice(100,225) ,\
               "B0_0.1" :slice(80,210)  ,\
               }
    # runFourierModesSliced(sliced=True)
    sub.setWalltime("00:10:00")
    # runGrowthRates()
    sub.setWalltime("00:05:00")
    # runEnergy(sliced=False)
    tSlices = {\
               "B0_0.02":None,\
               "B0_0.04":None,\
               "B0_0.06":slice(1200,None),\
               "B0_0.08":slice(1000,None),\
               "B0_0.1" :None,\
               }
    # runEnergy(sliced=True)
    sub.setWalltime("00:10:00")
    # runPosOfFluct()
    # runPSD2D()
    runSkewKurt()
