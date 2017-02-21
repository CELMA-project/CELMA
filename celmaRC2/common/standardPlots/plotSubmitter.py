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
from .analyticGrowthRates import analyticGrowthRatesPlot
from .blobs import (blobRadialFlux          ,\
                    blobWaitingTimePulsePlot,\
                    blobTimeTracesPlot      ,\
                    blob2DPlot)
from .blobDensPDF import blobDensPDF
from .combinedPlots import combinedPlotsPlot
from .fields1D import fields1DAnimation
from .fields2D import fields2DAnimation
from .fourierModes import fourierModesPlot
from .growthRates import growthRatesPlot
from .energy import energyPlot
from .growthRates import growthRatesPlot
from .magnitudeSpectrum import magnitudeSpectrumPlot
from .performance import performancePlot
from .phaseShift import phaseShiftPlot
from .posOfFluct import posOfFluctPlot
from .PSD2D import PSD2DPlot
from .skewKurt import skewKurtPlot
from .totalFlux import totalFluxPlot
from .zonalFlow import zonalFlowPlot
from copy import copy

#{{{PlotSubmitter
class PlotSubmitter(object):
    """Class used to submit the standard plots"""

    #{{{constructor
    def __init__(self, directory, scanParameter, boussinesq=False):
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
        boussinesq : bool
            Whether or not the boussinesq approximation is used
        """
        #}}}

        self._scanParameter = scanParameter

        # Set the folders to use
        with open(os.path.join(directory, "dmpFoldersDict.pickle"), "rb") as f:
                self._dmpFolders = pickle.load(f)
        allDir = ["init", "expand", "linear", "turbulence", "extraTurbulence"]
        fromLinear = ["linear", "turbulence", "extraTurbulence"]
        # Pop extra turbulence if not present
        if len(self._dmpFolders["extraTurbulence"]) == 0:
            allDir.remove("extraTurbulence")
            fromLinear.remove("extraTurbulence")
        allDir = tuple(allDir)
        fromLinear = tuple(fromLinear)
        self._mergeAll = pathMerger(self._dmpFolders, allDir)
        self._mergeInitAndExpand =\
           pathMerger(self._dmpFolders, ("init", "expand"))
        self._mergeFromLinear = pathMerger(self._dmpFolders, fromLinear)

        # Ranges to loop over
        self._paramKeys = tuple(sorted(list(self._mergeFromLinear.keys())))
        self._rangeJobs = range(len(self._paramKeys))

        # Generate the submitter
        self.sub = PBSSubmitter()
        self.sub.setNodes(nodes=1, ppn=20)
        self.sub.setQueue("xpresq")
        self.sub.setWalltime("00:15:00")

        # Create default plotSuperKwargs
        self._plotSuperKwargs = {\
                                 "showPlot"        : False        ,\
                                 "savePlot"        : True         ,\
                                 "savePath"        : None         ,\
                                 "savePathFunc"    : "onlyScan"   ,\
                                 "extension"       : None         ,\
                                 "timeStampFolder" : False        ,\
                                 # scanParameter needed in onlyScans
                                 "scanParameter"   : scanParameter,\
                                }

        # Set memeber data
        self._boussinesq = boussinesq
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

    #{{{updatePlotSuperKwargs
    def updatePlotSuperKwargs(self, updateDict):
        #{{{docstring
        """
        Function to update the plotSuperKwargs

        Paramteres
        ----------
        updateDict : dict
            Dictionary used to update the member plotSuperKwargs
        """
        #}}}

        self._plotSuperKwargs.update(updateDict)
    #}}}

    #{{{runBlobs
    def runBlobs(self         ,\
                 modes="all"  ,\
                 fluct="both" ,\
                 condition=3  ,\
                 varName ="n" ,\
                 phiCont=False,\
                 plotAll=False):
        #{{{docstring
        """
        Runs blobs.

        Parameters
        ----------
        modes : ["all"|"perp"|"par"|"pol"]
            2D mode to use
        fluct : ["both"|bool]
            Fluctuations to use
        condition : float
            The condition in the conditional average will be set to
            flux.std()*condition
        varName : str
            Variable to investigate
            NOTE: Condition is still on the radial density flux.
        phiCont : bool
            If True, phi contours will be overplotted on the 2D plots.
        plotAll : bool
           If True: The individual blobs will be plotted.
        """
        #}}}

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)

        if modes == "all":
            modes = ("perp", "par", "pol")
        else:
            modes = (modes,)
        if fluct == "both":
            flucts = (True, False)
        else:
            flucts = (fluct,)

        for dmp_folders, key, nr in loopOver:

            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders   = (dmp_folders,)

            args = (dmp_folders          ,\
                    collectPaths         ,\
                    self._plotSuperKwargs,\
                    tSlice               ,\
                    varName              ,\
                    condition            ,\
                    phiCont              ,\
                    plotAll              ,\
                    )

            kwargs = {}
            self.sub.setJobName("blobRadialFlux{}".format(nr))
            self.sub.submitFunction(blobRadialFlux,\
                                    args=args, kwargs=kwargs)

            self.sub.setJobName("blobWaitingTimePulse{}".format(nr))
            self.sub.submitFunction(blobWaitingTimePulsePlot,\
                                    args=args, kwargs=kwargs)

            self.sub.setJobName("blobTimeTrace{}".format(nr))
            self.sub.submitFunction(blobTimeTracesPlot,\
                                    args=args, kwargs=kwargs)

            for mode in modes:
                for b in flucts:
                    kwargs = {"mode":mode, "fluct":b}
                    if b:
                        fluct = "-fluct"
                    else:
                        fluct = ""
                    self.sub.setJobName("blob2DPlot-{}{}-{}".\
                                        format(mode,fluct,nr))
                    self.sub.submitFunction(blob2DPlot,\
                                            args=args, kwargs=kwargs)
    #}}}

    #{{{runBlobDensPDF
    def runBlobDensPDF(self):
        """
        Runs the density PDFs for the blobs
        """
        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("blobDensPDF{}".format(nr))
            self.sub.submitFunction(blobDensPDF, args=args, kwargs=kwargs)
    #}}}

    #{{{runCominedPlots
    def runCominedPlots(self):
        """
        Runs the combined plots
        """
        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("combinedPlotsSliced{}".format(nr))
            self.sub.submitFunction(combinedPlotsPlot,args=args,kwargs=kwargs)
    #}}}

    #{{{runAnalyticGrowthRates
    def runAnalyticGrowthRates(self):
        """
        Runs the growth rates
        """

        # NOTE: The ordering of param is in descending order (because of the
        #       organization in PBSScan)
        dmp_folders = (self._mergeFromLinear["param0"][-1],)
        steadyStatePaths = self._dmpFolders["expand"]

        # Local modification of plotSuperKwargs
        plotSuperKwargs = copy(self._plotSuperKwargs)
        newVals = {"savePath" : "all", "savePathFunc" : None}
        plotSuperKwargs.update(newVals)

        args = (dmp_folders        ,\
                steadyStatePaths   ,\
                self._scanParameter,\
                plotSuperKwargs    ,\
                )
        self.sub.setJobName("analyticGrowthRates")
        self.sub.submitFunction(analyticGrowthRatesPlot, args=args)
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
        loopOver = zip(self._dmpFolders["turbulence"],\
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
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            if sliced:
                kwargs = {"tSlice":tSlice}
                self.sub.setJobName("energySliced{}".format(nr))
            else:
                kwargs = {}
                self.sub.setJobName("energy{}".format(nr))
            self.sub.submitFunction(energyPlot, args=args, kwargs=kwargs)
    #}}}

    #{{{runFields1DAnim
    def runFields1DAnim(self, hyperIncluded=False):
        #{{{docstring
        """
        Runs the fields 1D for init and expand

        Parameters
        ----------
        hyperIncluded : bool
            If hyper viscosities are used.
        """
        #}}}

        loopOver = zip(self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            collectPaths = self._mergeInitAndExpand[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            kwargs = {"hyperIncluded" : hyperIncluded,\
                      "boussinesq"    : self._boussinesq}
            self.sub.setJobName("fields1D{}".format(nr))
            self.sub.submitFunction(fields1DAnimation, args=args, kwargs=kwargs)
    #}}}

    #{{{runFields2DAnim
    def runFields2DAnim(self, varName = "n", fluct=False):
        #{{{docstring
        """
        Runs the fields 2D for linear and turbulence

        Parameters
        ----------
        varName : str
            Variable to animate.
        fluct : bool
            Whether or not to plot the fluctuations.
        """
        #}}}

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)
            kwargs = {"varName":varName, "fluct":fluct}
            if fluct:
                self.sub.setJobName("fields2Dfluct{}".format(nr))
            else:
                self.sub.setJobName("fields2D{}".format(nr))
            self.sub.submitFunction(fields2DAnimation, args=args, kwargs=kwargs)
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

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            if sliced:
                tSlice = self._findSlices(dmp_folders, self._linearTSlices)
                if tSlice is None:
                    continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders   = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)
            if sliced:
                kwargs = {"tSlice":(tSlice,)}
                self.sub.setJobName("fourierModesSliced{}".format(nr))
            else:
                kwargs = {}
                self.sub.setJobName("fourierModes{}".format(nr))
            self.sub.submitFunction(fourierModesPlot, args=args, kwargs=kwargs)
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
        steadyStatePaths = self._dmpFolders["expand"]

        # NOTE: The ordering is of the keys are in descending order (because
        #       of the organization in PBSScan)
        growthTSlicesKeys =\
                tuple(sorted(list(self._linearTSlices.keys()), reverse=True))
        growthTSlices =\
                tuple(self._linearTSlices[key] for key in growthTSlicesKeys)

        # Local modification of plotSuperKwargs
        plotSuperKwargs = copy(self._plotSuperKwargs)
        newVals = {"savePath" : "all", "savePathFunc" : None}
        plotSuperKwargs.update(newVals)

        args = (dmp_folders        ,\
                scanCollectPaths   ,\
                steadyStatePaths   ,\
                self._scanParameter,\
                growthTSlices      ,\
                plotSuperKwargs    ,\
                )
        self.sub.setJobName("growthRates")
        self.sub.submitFunction(growthRatesPlot, args=args)
    #}}}

    #{{{runMagnitudeSpectrum
    def runMagnitudeSpectrum(self):
        """
        Runs the magnitude spectrum
        """

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            tSlice = self._findSlices(dmp_folders, self._linearTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders   = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)

            kwargs = {}
            self.sub.setJobName("magnitudeSpectrum{}".format(nr))
            self.sub.submitFunction(magnitudeSpectrumPlot, args=args, kwargs=kwargs)
    #}}}

    #{{{runSnapShotsSameScanVal
    def runSnapShotsSameScanVal(self, paramKey, slices,\
                                fluct = False, varName = "n", vMaxVMin = None,\
                                yInd = 16):
        #{{{docstring
        """
        Gets snapshot of a specific scan value.

        Parameters
        ----------
        paramKey : str
            What paramKey to use collective collect from
        slices : tuple
            Tuple of the tSlices to use.
        fluct : bool
            Whether or not to display the fluctuations.
        varName : str
            Variable to animate.
        vMaxVMin : [None|tuple]
            Max and min of the colors (if not None)
        yInd : int
            Parallel index to slice at.
        """
        #}}}

        plotSuperKwargs = copy(self._plotSuperKwargs)

        if vMaxVMin is not None:
            plotSuperKwargs["vmax"] = (vMaxVMin[0],)
            plotSuperKwargs["vmin"] = (vMaxVMin[1],)
            plotSuperKwargs["levels"] = None

        collectPaths    = self._mergeFromLinear[paramKey]
        steadyStatePath = self._dmpFolders["expand"][0]
        dmp_folders  = (collectPaths[-1],)
        args = (dmp_folders,\
                collectPaths,\
                steadyStatePath,\
                plotSuperKwargs)
        for nr, tSlice in enumerate(slices):
            kwargs = {"varName":varName,\
                      "fluct":fluct,\
                      "tSlice":tSlice,
                      "yInd" : yInd}
            self.sub.setJobName("snapShotsSameScanVal{}".format(nr))
            self.sub.submitFunction(fields2DAnimation, args=args, kwargs=kwargs)
            # Sleep to ensure that tmp files will have different names
            # FIXME: This is no guaranty for different names, but a
            #        workaround the cascade of variables arguments tru
            #        the system.
            #        Will fail if there is queue on the system
            sleep(30)
    #}}}

    #{{{runSnapShotDifferentScanVals
    def runSnapShotDifferentScanVals(self, slices,\
                                     fluct = True, varName = "n", yInd = 16):
        #{{{docstring
        """
        Gets snapshot of all the scan values.

        Parameters
        ----------
        slices : dict
            Dictionary containing the tSlices.
        fluct : bool
            Whether or not to display the fluctuations.
        varName : str
            Variable to animate.
        yInd : str
            Parallel index to slice at.
        """
        #}}}

        fluct = True

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            tSlice = self._findSlices(dmp_folders, slices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)
            kwargs = {"varName":varName,\
                      "fluct":fluct,\
                      "tSlice":tSlice,
                      "yInd" : yInd}
            self.sub.setJobName("snapShotDifferentScanVals{}".format(nr))
            self.sub.submitFunction(fields2DAnimation, args=args, kwargs=kwargs)
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
        loopOver = zip(self._dmpFolders["turbulence"],\
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
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            self.sub.submitFunction(performancePlot, args=args, kwargs=kwargs)
    #}}}

    #{{{runPhaseShift
    def runPhaseShift(self):
        """
        Runs the phase shift
        """

        # NOTE: The ordering of param is in descending order (because of the
        #       organization in PBSScan)
        dmp_folders = (self._mergeFromLinear["param0"][-1],)
        keys = tuple(sorted(list(self._mergeFromLinear.keys())))
        scanCollectPaths = tuple(self._mergeFromLinear[key] for key in keys)
        steadyStatePaths = self._dmpFolders["expand"]

        # NOTE: The ordering is of the keys are in descending order (because
        #       of the organization in PBSScan)
        growthTSlicesKeys =\
                tuple(sorted(list(self._linearTSlices.keys()), reverse=True))
        growthTSlices =\
                tuple(self._linearTSlices[key] for key in growthTSlicesKeys)

        # Local modification of plotSuperKwargs
        plotSuperKwargs = copy(self._plotSuperKwargs)
        newVals = {"savePath" : "all", "savePathFunc" : None}
        plotSuperKwargs.update(newVals)

        args = (dmp_folders        ,\
                scanCollectPaths   ,\
                steadyStatePaths   ,\
                self._scanParameter,\
                growthTSlices      ,\
                plotSuperKwargs    ,\
                )
        self.sub.setJobName("phaseShift")
        self.sub.submitFunction(phaseShiftPlot, args=args)
    #}}}

    #{{{runPosOfFluct
    def runPosOfFluct(self):
        """
        Runs the position of fluct
        """

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("posOfFluctSliced{}".format(nr))
            self.sub.submitFunction(posOfFluctPlot, args=args, kwargs=kwargs)
    #}}}

    #{{{runPSD2D
    def runPSD2D(self):
        """
        Runs the PSD2D
        """

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("PSD2DPlotSliced{}".format(nr))
            self.sub.submitFunction(PSD2DPlot, args=args, kwargs=kwargs)
    #}}}

    #{{{runSkewKurt
    def runSkewKurt(self):
        """
        Runs the skewness and kurtosis
        """

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("skewnessKurtosisSliced{}".format(nr))
            self.sub.submitFunction(skewKurtPlot, args=args, kwargs=kwargs)
    #}}}

    #{{{runSteadyState
    def runSteadyState(self):
        """
        Runs the 1D steady state profiles
        """

        hyperIncluded = False
        tSlice = slice(-1,-1)

        loopOver = zip(self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:

            collectPaths = self._mergeInitAndExpand[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            kwargs = {"hyperIncluded" : hyperIncluded,\
                      "boussinesq"    : self._boussinesq,\
                      "tSlice"        : tSlice,\
                     }
            self.sub.setJobName("steadyState{}".format(nr))
            self.sub.submitFunction(fields1DAnimation, args=args, kwargs=kwargs)
    #}}}

    #{{{runTotalFlux
    def runTotalFlux(self):
        """
        Runs the total flux
        """

        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, key, nr in loopOver:
            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders, collectPaths, self._plotSuperKwargs)
            self.sub.setJobName("totalFlux{}".format(nr))
            self.sub.submitFunction(totalFluxPlot, args=args)
    #}}}

    #{{{runZonalFlow
    def runZonalFlow(self):
        """
        Runs the zonal flow
        """
        loopOver = zip(self._dmpFolders["turbulence"],\
                       self._dmpFolders["expand"],\
                       self._paramKeys,\
                       self._rangeJobs)
        for dmp_folders, steadyStatePath, key, nr in loopOver:

            # Find tSlice
            tSlice = self._findSlices(dmp_folders, self._satTurbTSlices)
            if tSlice is None:
                continue

            collectPaths = self._mergeFromLinear[key]
            dmp_folders  = (dmp_folders,)
            args = (dmp_folders,\
                    collectPaths,\
                    steadyStatePath,\
                    self._plotSuperKwargs)
            kwargs = {"tSlice":tSlice}
            self.sub.setJobName("zonalFlowSliced{}".format(nr))
            self.sub.submitFunction(zonalFlowPlot, args=args, kwargs=kwargs)
    #}}}
#}}}
