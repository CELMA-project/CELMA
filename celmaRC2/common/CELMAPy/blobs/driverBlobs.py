#!/usr/bin/env python

"""
Contains drivers for the blobs
"""

from ..plotHelpers import getVmaxVminLevels
from ..superClasses import DriverSuperClass
from ..fields2D import (PlotAnim2DPerp,\
                        PlotAnim2DPar,\
                        PlotAnim2DPol,\
                       )
from .collectAndCalcBlobs import CollectAndCalcBlobs
from .plotBlobs import PlotTemporalStats, PlotBlobTimeTrace
from multiprocessing import Process

#{{{prepareBlobs
def prepareBlobs(collectPaths     ,\
                 slices           ,\
                 pctPadding       ,\
                 convertToPhysical,\
                ):
    #{{{docstring
    """
    Driver for plotting blobs.

    Parameters
    ----------
    collectPaths : tuple
        Tuple from where to collect
    slices : tuple of tuples
        Tuple the indices to use.
        On the form (xInd, yInd, zInd, tSlice)
    pctPadding : float
        Padding around the maximum pulsewidth which satisfies the
        condition.
        Measured in percent.
    convertToPhysical : bool
        Whether or not to convert to physical

    Returns
    -------
    ccb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    """
    #}}}

    ccb = CollectAndCalcBlobs(collectPaths, slices, convertToPhysical)

    ccb.prepareCollectAndCalc()

    return ccb
#}}}

#{{{driverRadialFluxWStdLines
def driverRadialFluxWStdLines(ccb            ,\
                              plotSuperKwargs,\
                              stdConditions  ,\
                             ):
    #{{{docstring
    """
    Driver for plotting radial fluxes with standard deviation lines.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    stdConditions : tuple
        Conditions for the conditional averages.
    """
    #}}}

    radialFlux = getRadialFlux()

    # Plot
    ptt = PlotRadialFlux(ccb.uc, **plotSuperKwargs)
    ptt.setData(radialFlux, mode, stdLines=stdConditions)
    ptt.plotSaveShowRadialFlux()
#}}}

#{{{driverWaitingTimePulse
def driverWaitingTimePulse(ccb, plotSuperKwargs, normed = False):
    #{{{docstring
    """
    Driver which plots the waiting time and pulse width statistics.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    normed : bool
        Wheter or not the histogram should be normed
    """
    #}}}

    holesWaitingTime, holesPulseWidths =\
        ccb.getWaitingTimesAndPulseWidth("holes")
    blobsWaitingTime, blobsPulseWidths =\
        ccb.getWaitingTimesAndPulseWidth("blobs")

    pts = PlotTemporalStats(ccb.uc          ,\
                            **plotSuperKwargs)

    if len(blobsPulseWidths) > 0:
        pts.setData(blobsWaitingTime, blobsPulseWidths, "blobs", normed)
        pts.plotSaveShowTemporalStats()
    else:
        print("No blob time statistic made as no blobs were detected")

    if len(holesPulseWidths) > 0:
        pts.setData(holesWaitingTime, holesPulseWidths, "holes", normed)
        pts.plotSaveShowTemporalStats()
    else:
        print("No hole time statistic made as no holes were detected")
#}}}

#{{{driverBlobTimeTraces
def driverBlobTimeTraces(ccb, plotSuperKwargs):
    #{{{docstring
    """
    Driver which plots the time traces.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    timeTraceBlobsAvg, timeTraceBlobs, timeTraceHolesAvg, timeTraceHoles =\
        ccb.executeCollectAndCalc1D()

    pbtt = PlotBlobTimeTrace(ccb.uc, **plotSuperKwargs)

    # Blobs
    if len(timeTraceBlobs) > 0:
        for theDict in (timeTraceBlobsAvg, *timeTraceBlobs):
            pbtt.setData(theDict, "blobs")
            pbtt.plotSaveShowTimeTrace()
    else:
        print("No blobs plotted as no blobs were detected")

    # Holes
    if len(timeTraceHoles) > 0:
        for theDict in (timeTraceHolesAvg, *timeTraceHoles):
            pbtt.setData(theDict, "holes")
            pbtt.plotSaveShowTimeTrace()
    else:
        print("No holes plotted as no holes were detected")
#}}}

#{{{get2DData
def get2DData(ccb, mode, fluct):
    #{{{docstring
    """
    Driver which collects the 2D slices.

    Exists as an own entity if combined plots should be made.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    mode : ["perp"|"par"|"pol"]
        The mode to collect.
    fluct : bool
        Whether or not the fluctuations will be collected.

    Returns
    -------
    blobs2DAvg : dict
        Dictionary of the averaged blob.
        Contains the keys:
            * "n"    - The 2D variable.
            * "nPPi" - The 2D variable pi away from the set zInd
                       (only when mode is "par").
            * "time" - The corresponding time.
            * "X"    - The X-mesh.
            * "Y"    - The Y-mesh.
            * pos    - The position of the fixed index
    blobs2D : tuple
        Tuple containing the dictionaries used to calculate the
        averaged blob.
        Each element contains the same keys as blobs2DAvg.
    holes2DAvg : dict
        Dictionary of the averaged hole.
        Contains the same keys as blobs2DAvg.
    holes2D : tuple
        Tuple containing the dictionaries used to calculate the
        averaged blob.
        Each element contains the same keys as blobs2DAvg.
    """
    #}}}

    blobs2DAvg, blobs2D, holes2DAvg, holes2D =\
        ccb.executeCollectAndCalc2D(mode, fluct)

    return blobs2DAvg, blobs2D, holes2DAvg, holes2D
#}}}

#{{{driverPlot2DData
def driverPlot2DData(ccb, mode, fluct, plotSuperKwargs):
    #{{{docstring
    """
    Driver which plots the 2D data.

    NOTE: mode == "par" is slower then the rest, as this requires more
          opening of files.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    mode : ["perp"|"par"|"pol"]
        The mode to collect.
    fluct : bool
        Whether or not the fluctuations will be collected.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    varName = "n"
    # We would like all the frames to have the same max/min
    varyMaxMin = False

    blobs2DAvg, blobs2D, holes2DAvg, holes2D =\
        get2DData(ccb, mode, fluct)

    blobsAndHoles = ((blobs2DAvg, *blobs2D), (holes2DAvg, *holes2D))
    for curTuple, blobOrHole in zip(blobsAndHoles,("blobs", "holes")):
        if len(curTuple[0]) == 0:
            continue
        plotSuperKwargs["blobOrHole"] = blobOrHole
        for nr, blob in enumerate(curTuple):
            if nr == 0:
                plotSuperKwargs["averagedBlobOrHole"] = True
            else:
                plotSuperKwargs["averagedBlobOrHole"] = False
            args = (blob          ,\
                    varName       ,\
                    varyMaxMin    ,\
                    fluct         ,\
                    ccb           ,\
                    plotSuperKwargs)
            if mode == "perp":
                plotBlob2DPerp(*args)
            elif mode == "par":
                plotBlob2DPar(*args)
            elif mode == "pol":
                plotBlob2DPol(*args)
#}}}

#{{{plotBlob2DPerp
def plotBlob2DPerp(blob, varName, varyMaxMin, fluct, ccb, plotSuperKwargs):
    #{{{docstring
    """
    Performs the plotting of the blob in a 2D perp plane.

    Parameters
    -----------
    blob : dict
        Dictionary with the variables to be plotted.
        Contains the keys:
            * "n"    - The 2D variable.
            * "nPPi" - The 2D variable pi away from the set zInd
                       (only when mode is "par").
            * "time" - The corresponding time.
            * "X"    - The X-mesh.
            * "Y"    - The Y-mesh.
            * pos    - The position of the fixed index
    varName : str
        Name of the variable.
    varyMaxMin : bool
        Whether or not to vary the max and min for each plot in the
        colorbar.
    fluct : bool
        Wheter or not only the fluctuations are given as an input
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the PlotSuperClass class.
    """
    #}}}
    # Set the plot limits
    tupleOfArrays = (blob[varName],)
    vmax, vmin, levels =\
            getVmaxVminLevels(plotSuperKwargs,\
                              tupleOfArrays,\
                              fluct,\
                              varyMaxMin)

    # Create the plotting object
    p2DPerp = PlotAnim2DPerp(ccb.uc          ,\
                             fluct    = fluct,\
                             **plotSuperKwargs)
    p2DPerp.setContourfArguments(vmax, vmin, levels)
    p2DPerp.setPerpData(blob["X"],\
                        blob["Y"],\
                        blob[varName],\
                        blob["time"],\
                        blob["zPos"],\
                        varName)
    p2DPerp.plotAndSavePerpPlane()
#}}}

#{{{plotBlob2DPar
def plotBlob2DPar(blob, varName, varyMaxMin, fluct, ccb, plotSuperKwargs):
    #{{{docstring
    """
    Performs the plotting of the blob in a 2D par plane.

    Parameters
    -----------
    blob : dict
        Dictionary with the variables to be plotted.
        Contains the keys:
            * "n"    - The 2D variable.
            * "nPPi" - The 2D variable pi away from the set zInd
                       (only when mode is "par").
            * "time" - The corresponding time.
            * "X"    - The X-mesh.
            * "Y"    - The Y-mesh.
            * pos    - The position of the fixed index
    varName : str
        Name of the variable.
    varyMaxMin : bool
        Whether or not to vary the max and min for each plot in the
        colorbar.
    fluct : bool
        Wheter or not only the fluctuations are given as an input
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the PlotSuperClass class.
    """
    #}}}
    # Set the plot limits
    tupleOfArrays = (blob[varName], blob[varName+"PPi"])
    vmax, vmin, levels =\
            getVmaxVminLevels(plotSuperKwargs,\
                              tupleOfArrays,\
                              fluct,\
                              varyMaxMin)

    # Create the plotting object
    p2DPar = PlotAnim2DPar(ccb.uc          ,\
                           fluct    = fluct,\
                           **plotSuperKwargs)
    p2DPar.setContourfArguments(vmax, vmin, levels)
    p2DPar.setParData(blob["X"]          ,\
                      blob["Y"]          ,\
                      blob[varName]      ,\
                      blob[varName+"PPi"],\
                      blob["time"]       ,\
                      blob["thetaPos"]   ,\
                      varName)
    p2DPar.plotAndSaveParPlane()
#}}}

#{{{plotBlob2DPol
def plotBlob2DPol(blob, varName, varyMaxMin, fluct, ccb, plotSuperKwargs):
    #{{{docstring
    """
    Performs the plotting of the blob in a 2D pol plane.

    Polameters
    -----------
    blob : dict
        Dictionary with the variables to be plotted.
        Contains the keys:
            * "n"    - The 2D variable.
            * "nPPi" - The 2D variable pi away from the set zInd
                       (only when mode is "pol").
            * "time" - The corresponding time.
            * "X"    - The X-mesh.
            * "Y"    - The Y-mesh.
            * pos    - The position of the fixed index
    varName : str
        Name of the variable.
    varyMaxMin : bool
        Whether or not to vary the max and min for each plot in the
        colorbar.
    fluct : bool
        Wheter or not only the fluctuations are given as an input
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the PlotSuperClass class.
    """
    #}}}
    # Set the plot limits
    tupleOfArrays = (blob[varName],)
    vmax, vmin, levels =\
            getVmaxVminLevels(plotSuperKwargs,\
                              tupleOfArrays,\
                              fluct,\
                              varyMaxMin)

    # Create the plotting object
    p2DPol = PlotAnim2DPol(ccb.uc          ,\
                             fluct    = fluct,\
                             **plotSuperKwargs)
    p2DPol.setContourfArguments(vmax, vmin, levels)
    p2DPol.setPolData(blob["X"],\
                      blob["Y"],\
                      blob[varName],\
                      blob["time"],\
                      blob["rhoPos"],\
                      varName)
    p2DPol.plotAndSavePolPlane()
#}}}

#{{{DriverBlobs
class DriverBlobs(DriverSuperClass):
    """
    Class for driving of the plotting of the blobs.
    """
    #{{{Constructor
    def __init__(self                       ,\
                 dmp_folders                ,\
                 slices                     ,\
                 pctPadding                 ,\
                 convertToPhysical          ,\
                 plotSuperKwargs            ,\
                 stdConditions = (4,3,2.5,2),\
                 normed        = False      ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Set the member data
            * Updates the plotSuperKwargs

        Parameters
        ----------
        dmp_folders : tuple
            Tuple of the dmp_folder (output from bout_runners).
        slices : tuple of tuples
            Tuple the indices to use.
            On the form (xInd, yInd, zInd, tSlice)
        pctPadding : float
            Padding around the maximum pulsewidth which satisfies the
            condition.
            Measured in percent.
        convertToPhysical : bool
            Whether or not to convert to physical
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        stdConditions : tuple
            Conditions for the conditional averages.
        normed : bool
            Wheter or not to norm the histogram
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._slices            = slices
        self._convertToPhysical = convertToPhysical
        self._pctPadding        = pctPadding
        self._stdConditions     = stdConditions
        self._normed            = normed

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"blobs"})
        self._plotSuperKwargs = plotSuperKwargs

        # Prepare the blobs
        self._ccb =\
            prepareBlobs(self._collectPaths,\
                         slices            ,\
                         pctPadding        ,\
                         convertToPhysical)
    #}}}

    #{{{driverAll
    def driverAll(self):
        #{{{docstring
        """
        Runs all the drivers
        """
        #}}}

        self.driverWaitingTimePulse()
        self.driverBlobTimeTraces()
        self.driverRadialFluxWStdLines()
        modes = ("perp", "par", "pol")
        for mode in modes:
            self.setMode(mode)
            for b in (True, False):
                self.setFluct(b)
                self.driverPlot2DData()
    #}}}

    #{{{driverRadialFluxWStdLines
    def driverRadialFluxWStdLines(self):
        #{{{docstring
        """
        Wrapper to driverRadialFluxWStdLines
        """
        #}}}

        args = (\
                ccb            ,\
                plotSuperKwargs,\
                stdConditions  ,\
               )
        if self._useSubProcess:
            processes = Process(target = driverRadialFluxWStdLines,\
                                args   = args                     )
            processes.start()
        else:
            driverRadialFluxWStdLines(*args, **kwargs)
    #}}}

    #{{{driverWaitingTimePulse
    def driverWaitingTimePulse(self):
        #{{{docstring
        """
        Wrapper to driverWaitingTimePulse
        """
        #}}}

        args   = (self._ccb, self._plotSuperKwargs)
        kwargs = {"normed":self._normed}
        if self._useSubProcess:
            processes = Process(target = driverWaitingTimePulse,\
                                args   = args                  ,\
                                kwargs = kwargs)
            processes.start()
        else:
            driverWaitingTimePulse(*args, **kwargs)
    #}}}

    #{{{driverBlobTimeTraces
    def driverBlobTimeTraces(self):
        #{{{docstring
        """
        Wrapper to driverBlobTimeTraces
        """
        #}}}

        args  = (self._ccb, self._plotSuperKwargs)
        if self._useSubProcess:
            processes = Process(target = driverBlobTimeTraces,\
                                args   = args                ,\
                               )
            processes.start()
        else:
            driverBlobTimeTraces(*args)
    #}}}

    #{{{driverPlot2DData
    def driverPlot2DData(self):
        #{{{docstring
        """
        Wrapper to driverPlot2DData
        """
        #}}}

        args  = (self._ccb            ,\
                 self._mode           ,\
                 self._fluct          ,\
                 self._plotSuperKwargs,\
                )
        if self._useSubProcess:
            processes = Process(target = driverPlot2DData,\
                                args   = args            ,\
                               )
            processes.start()
        else:
            driverPlot2DData(*args)
    #}}}

    #{{{setMode
    def setMode(self, mode):
        #{{{docstring
        """
        Setter for self._mode
        """
        #}}}
        implemeted = ("perp" ,"par", "pol")
        if not(mode in implemeted):
            message = "Got '{}', expected one of the following:\n".format(mode)
            for i in implemented:
                message += "{}".format(i)
            raise ValueError(message)

        self._mode = mode
    #}}}

    #{{{setFluct
    def setFluct(self, fluct):
        #{{{docstring
        """
        Setter for self._fluct
        """
        #}}}
        self._fluct = fluct
    #}}}
##}}}
