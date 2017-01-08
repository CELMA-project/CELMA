#!/usr/bin/env python

"""Contains the PlotSubmitter class."""

import pickle
from time import sleep

import os, sys
# If we add to sys.path, then it must be an absolute path
commonDir = os.path.abspath("./../")
# Sys path is a list of system paths
sys.path.append(commonDir)

from CELMAPy.driverHelpers import PBSSubmitter, pathMerger
from .combinedPlots import combinedPlotsPlot
from .fourierModes import fourierModesPlot
from .growthRates import growthRatesPlot
from .energy import energyPlot
from .performance import performancePlot
from .posOfFluct import posOfFluctPlot
from .PSD2D import PSD2DPlot
from .skewKurt import skewKurtPlot
from .zonalFlow import zonalFlowPlot

#{{{PlotSubmitter
class PlotSubmitter(object):
    """Class used to submit the standard plots"""

    #{{{constructor
    def __init__(self, directory, scanParameter):
        #{{{docstring
        """
        Constructor for the PlotSubmitter class.

        The constructor will:
            * Load the pickle containing the dmp_folders.
            * Make variables to loop over
            * Generate the submitter to be used.

        Parameters
        ----------
        directory : str
            Root directory of the simulations.
        scanParameter : str
            The scan parameter.
        """
        #}}}

        self._scanParameter = scanParameter

        with open(os.path.join(directory, "dmpFoldersDict.pickle"), "rb") as f:
                dmpFolders = pickle.load(f)

        mergeFromLinear =\
            pathMerger(dmpFolders, ("linear", "turbulence", "extraTurbulence"))

        mergeAll =\
           pathMerger(dmpFolders,\
                ("init", "expand", "linear", "turbulence", "extraTurbulence"))

        paramKeys = tuple(sorted(list(mergeFromLinear.keys())))
        rJobs     = range(len(paramKeys))

        # Generate the submitter
        sub = PBSSubmitter()
        sub.setNodes(nodes=1, ppn=2)
        sub.setQueue("xpresq")
        sub.setWalltime("00:15:00")
    #}}}

    #{{{_findSlices
    def _findSlices(self, dmp_folders):
        #{{{docstring
        """
        Finds the temporal slices to use for the dmp_folder.

        Parameters
        ----------
        dmp_folders : str
            Dump folder under investigation

        Returns
        -------
        tSlice : slice
            The tSlice corresponding to the dmp_folder.
        """
        #}}}

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

    #{{{setLinearPhaseTSlices
    def setLinearPhaseTSlices(self, tSlices):
        """
        Set the time slices for the linear phase
        """

        self._linearTSlices = tSlices
    #}}}

    #{{{setSatTurbTSlices
    def setSatTurbTSlices(self, tSlices):
        """
        Set the time slices for the saturated turbulence phase
        """

        self._satTurbTSlices = tSlices
    #}}}

    #{{{runFourierModes
    def runFourierModes(self, sliced = False):
        #{{{docstring
        """
        Runs the fourier modes

        Parameters
        ----------
        sliced : bool
            Whether or not to slice the time
        """
        #}}}

        loopOver = zip(dmpFolders["turbulence"],\
                       dmpFolders["expand"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            if sliced:
                tSlice = _findSlices(dmp_folders)
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
    def runGrowthRates(self):
        """
        Runs the growth rates
        """

        # NOTE: The ordering of param is in descending order (because of the
        #       organization in PBSScan)
        dmp_folders = (mergeFromLinear["param0"][-1],)
        keys = tuple(sorted(list(mergeFromLinear.keys())))
        scanCollectPaths = tuple(mergeFromLinear[key] for key in keys)
        steadyStatePaths = dmpFolders["expand"]

        # NOTE: The ordering is of the keys are in descending order (because
        #       of the organization in PBSScan)
        growthTSlicesKeys = tuple(sorted(list(tSlices.keys()), reverse=True))
        growthTSlices = tuple(tSlices[key] for key in growthTSlicesKeys)

        args = (dmp_folders     ,\
                scanCollectPaths,\
                steadyStatePaths,\
                scanParameter   ,\
                growthTSlices)
        sub.setJobName("growthRates")
        sub.submitFunction(growthRatesPlot, args=args)
    #}}}

    #{{{runEnergy
    def runEnergy(self, sliced=False):
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
                tSlice = _findSlices(dmp_folders)
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
    def runPosOfFluct(self):
        """
        Runs the position of fluct
        """

        loopOver = zip(dmpFolders["turbulence"],\
                       dmpFolders["expand"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = _findSlices(dmp_folders)
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
    def runPSD2D(self):
        """
        Runs the PSD2D
        """

        loopOver = zip(dmpFolders["turbulence"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = _findSlices(dmp_folders)
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
    def runSkewKurt(self):
        """
        Runs the skewness and kurtosis
        """

        loopOver = zip(dmpFolders["turbulence"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = _findSlices(dmp_folders)
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

    #{{{runZonalFlow
    def runZonalFlow(self):
        """
        Runs the skewness and kurtosis
        """
        loopOver = zip(dmpFolders["turbulence"],\
                       dmpFolders["expand"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = _findSlices(dmp_folders)
            if tSlice is None:
                continue

            collectPaths = mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, steadyStatePath)
            kwargs = {"tSlice":tSlice}
            sub.setJobName("zonalFlowSliced{}".format(nr))
            sub.submitFunction(zonalFlowPlot, args=args, kwargs=kwargs)
            # Sleep 1.5 seconds to ensure that tmp files will have different names
            sleep(1.5)
    #}}}

    #{{{runCominedPlots
    def runCominedPlots(self):
        """
        Runs the combined plots
        """
        loopOver = zip(dmpFolders["turbulence"],\
                       dmpFolders["expand"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = _findSlices(dmp_folders)
            if tSlice is None:
                continue

            collectPaths = mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, steadyStatePath)
            kwargs = {"tSlice":tSlice}
            sub.setJobName("combinedPlotsSliced{}".format(nr))
            sub.submitFunction(combinedPlotsPlot, args=args, kwargs=kwargs)
            # Sleep 1.5 seconds to ensure that tmp files will have different names
            sleep(1.5)
    #}}}

    #{{{runPerformance
    def runPerformance(self, allFolders=False):
        #{{{docstring
        """
        Runs the performance plots

        Parameters
        ----------
        allFolders : bool
            If "init" and "expand" should be included in the plot.
        """
        #}}}
        loopOver = zip(dmpFolders["turbulence"],\
                       paramKeys,\
                       rJobs)
        for dmp_folders, key, nr in loopOver:

            if allFolders:
                collectPaths = mergeAll[key]
                sub.setJobName("performanceAll{}".format(nr))
                kwargs = {"allFolders":True}
            else:
                collectPaths = mergeFromLinear[key]
                sub.setJobName("performance{}".format(nr))
                kwargs = {"allFolders":False}
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths)
            sub.submitFunction(performancePlot, args=args, kwargs=kwargs)
            # Sleep 1.5 seconds to ensure that tmp files will have different names
            sleep(1.5)
    #}}}
#}}}
