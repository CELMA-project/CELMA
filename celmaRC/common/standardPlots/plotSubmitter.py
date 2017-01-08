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

        # Magic number
        self._sleepS = 1.5

        self._scanParameter = scanParameter

        with open(os.path.join(directory, "dmpFoldersDict.pickle"), "rb") as f:
                self._dmpFolderes = pickle.load(f)

        self._mergeFromLinear =\
            pathMerger(self._dmpFolderes,\
                       ("linear", "turbulence", "extraTurbulence"))

        self._mergeAll =\
           pathMerger(self._dmpFolderes,\
                ("init", "expand", "linear", "turbulence", "extraTurbulence"))

        self._paramKeys = tuple(sorted(list(self._mergeFromLinear.keys())))
        self._rangeJobs = range(len(self._paramKeys))

        # Generate the submitter
        self.sub = PBSSubmitter()
        self.sub.setNodes(nodes=1, ppn=2)
        self.sub.setQueue("xpresq")
        self.sub.setWalltime("00:15:00")
    #}}}

    #{{{_findSlices
    def _findSlices(self, dmp_folders, tSlices):
        #{{{docstring
        """
        Finds the temporal slices to use for the dmp_folder.

        Parameters
        ----------
        dmp_folders : str
            Dump folder under investigation
        tSlices : dict
            Dictionary containing the scanParameter

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

    #{{{runCominedPlots
    def runCominedPlots(self):
        """
        Runs the combined plots
        """
        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._dmpFolderes["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, steadyStatePath)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("combinedPlotsSliced{}".format(nr))
            self.sub.submitFunction(combinedPlotsPlot,args=args,kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
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
        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            if sliced:
                # Find tSlice
                tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
                if tSlice is None:
                    continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths)
            if sliced:
                kwargs = {"tSlice":tSlice}
                self.sub.setJobName("energySliced{}".format(nr))
            else:
                kwargs = {}
                self.sub.setJobName("energy{}".format(nr))
            self.sub.submitFunction(energyPlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
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

        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._dmpFolderes["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            if sliced:
                tSlice = self._findSlices(dmp_folders, self._linearTSlices)
                if tSlice is None:
                    continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders   = (dmp_folders,)
            args = (dmp_folders, collectPaths, steadyStatePath)
            if sliced:
                kwargs = {"tSlice":(tSlice,)}
                self.sub.setJobName("fourierModesSliced{}".format(nr))
            else:
                kwargs = {}
                self.sub.setJobName("fourierModes{}".format(nr))
            self.sub.submitFunction(fourierModesPlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
    #}}}

    #{{{runGrowthRates
    def runGrowthRates(self):
        """
        Runs the growth rates
        """

        # NOTE: The ordering of param is in descending order (because of the
        #       organization in PBSScan)
        dmp_folders = (self._mergeFromLinear["param0"][-1],)
        keys = tuple(sorted(list(self._mergeFromLinear.keys())))
        scanCollectPaths = tuple(self._mergeFromLinear[key] for key in keys)
        steadyStatePaths = self._dmpFolderes["expand"]

        # NOTE: The ordering is of the keys are in descending order (because
        #       of the organization in PBSScan)
        growthTSlicesKeys =\
                tuple(sorted(list(self._linearTSlices.keys()), reverse=True))
        growthTSlices =\
                tuple(self._linearTSlices[key] for key in growthTSlicesKeys)

        args = (dmp_folders        ,\
                scanCollectPaths   ,\
                steadyStatePaths   ,\
                self._scanParameter,\
                growthTSlices)
        self.sub.setJobName("growthRates")
        self.sub.submitFunction(growthRatesPlot, args=args)
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
        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            if allFolders:
                collectPaths = self._mergeAll[key]
                self.sub.setJobName("performanceAll{}".format(nr))
                kwargs = {"allFolders":True}
            else:
                collectPaths = self._mergeFromLinear[key]
                self.sub.setJobName("performance{}".format(nr))
                kwargs = {"allFolders":False}
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths)
            self.sub.submitFunction(performancePlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
    #}}}

    #{{{runPosOfFluct
    def runPosOfFluct(self):
        """
        Runs the position of fluct
        """

        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._dmpFolderes["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, steadyStatePath)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("posOfFluctSliced{}".format(nr))
            self.sub.submitFunction(posOfFluctPlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
    #}}}

    #{{{runPSD2D
    def runPSD2D(self):
        """
        Runs the PSD2D
        """

        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("PSD2DPlotSliced{}".format(nr))
            self.sub.submitFunction(PSD2DPlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
    #}}}

    #{{{runSkewKurt
    def runSkewKurt(self):
        """
        Runs the skewness and kurtosis
        """

        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("skewnessKurtosisSliced{}".format(nr))
            self.sub.submitFunction(skewKurtPlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
    #}}}

    #{{{runZonalFlow
    def runZonalFlow(self):
        """
        Runs the skewness and kurtosis
        """
        loopOver = zip(self._dmpFolderes["turbulence"],\
                       self._dmpFolderes["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, steadyStatePath)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("zonalFlowSliced{}".format(nr))
            self.sub.submitFunction(zonalFlowPlot, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            sleep(self._sleepS)
    #}}}
#}}}
