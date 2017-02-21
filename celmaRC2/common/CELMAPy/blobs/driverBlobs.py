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
from ..radialFlux import PlotRadialFlux
from ..PDF import CollectAndCalcPDF, PlotPDF
from ..superClasses import PlotSuperClass
from ..unitsConverter import UnitsConverter
from .collectAndCalcBlobs import CollectAndCalcBlobs
from .plotBlobs import (PlotTemporalStats,\
                        PlotBlobOrHoleTimeTraceSingle,\
                        PlotBlobAndHoleTimeTraceDouble)
from multiprocessing import Process
from collections import namedtuple
import os, pickle
import numpy as np

# This should be placed in the filescop in order for it to be pickable
Count = namedtuple("Count", ("blobs", "holes"))

#{{{prepareBlobs
def prepareBlobs(collectPaths     ,\
                 slices           ,\
                 pctPadding       ,\
                 convertToPhysical,\
                 condition = 3    ,\
                 picklePath = None,\
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
    condition : float
        The condition in the conditional average will be set to
        flux.std()*condition
    picklePath : [None|str]
        If set, the ccb will be pickled to the path if it doesn't
        exists, or read from the pickle if already exists

    Returns
    -------
    ccb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    """
    #}}}

    collect = True
    if picklePath:
        fileName = os.path.join(picklePath, "ccb.pickle")
        if os.path.exists(fileName):
            collect = False
            with open(fileName, "rb") as f:
                ccb = pickle.load(f)

    if collect:
        ccb = CollectAndCalcBlobs(collectPaths          ,\
                                  slices                ,\
                                  convertToPhysical     ,\
                                  condition =  condition,\
                                  )

        ccb.prepareCollectAndCalc()

    if picklePath and collect:
        with open(fileName, "wb") as f:
            pickle.dump(ccb, f, pickle.HIGHEST_PROTOCOL)

    return ccb
#}}}

#{{{driverRadialFlux
def driverRadialFlux(ccb            ,\
                     plotSuperKwargs,\
                    ):
    #{{{docstring
    """
    Driver which:

        1. Plots the time trace of the radial flux.
        2. Plots the PDF of the radial flux.
           NOTE: The labels will be wrong, but will be fixed in pickleTweaks
        3. Stores the blob count.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    radialFlux = ccb.getRadialFlux()

    # Plot the time trace
    ptt = PlotRadialFlux(ccb.uc, **plotSuperKwargs)
    ptt.setData(radialFlux, "fluct")
    ptt.plotSaveShowRadialFlux()

    # Plot the PDF
    # NOTE: The labels will be wrong, but will be fixed in pickleTweaks
    PDF = CollectAndCalcPDF.calcPDF(radialFlux)
    # Rename the keys for conversion dict to find the key
    for key in list(PDF.keys()):
        PDF[key]["fluxnPDFX"] = PDF[key].pop("nRadialFluxPDFX")
        PDF[key]["fluxnPDFY"] = PDF[key].pop("nRadialFluxPDFY")
    ppdf = PlotPDF(ccb.uc, **plotSuperKwargs)
    ppdf.setData(PDF, "fluct")
    ppdf.plotSaveShowPDF()

    # Store counts into a namedtuple
    blobCount, holeCount = ccb.getCounts()
    count = Count(blobCount, holeCount)
    # Pickle the named tuple
    fileName = os.path.join(ptt.getSavePath(), "counts.pickle")
    with open(fileName, "wb") as f:
        pickle.dump(count, f, pickle.HIGHEST_PROTOCOL)
        print("Counts saved to "+fileName)
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
def driverBlobTimeTraces(ccb, plotSuperKwargs, plotAll):
    #{{{docstring
    """
    Driver which plots the time traces.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    plotAll : bool
           If True: The individual blobs will be plotted.
    """
    #}}}

    timeTraceBlobsAvg, timeTraceBlobs, timeTraceHolesAvg, timeTraceHoles =\
        ccb.executeCollectAndCalc1D()

    if not(plotAll) and (len(timeTraceBlobs) > 0 and len(timeTraceHoles) > 0):
        # Plot the blobs and holes in one
        pbtt = PlotBlobAndHoleTimeTraceDouble(ccb.uc, **plotSuperKwargs)
        pbtt.setData(timeTraceBlobsAvg, timeTraceHolesAvg)
        pbtt.plotSaveShowTimeTrace()
    else:
        # Plot the blobs and holes separately
        pbtt = PlotBlobOrHoleTimeTraceSingle(ccb.uc, **plotSuperKwargs)
        # Blobs
        if len(timeTraceBlobs) > 0:
            if plotAll:
                loopOver = (timeTraceBlobsAvg, *timeTraceBlobs)
            else:
                loopOver = (timeTraceBlobsAvg, )
            for theDict in loopOver:
                pbtt.setData(theDict, "blobs")
                pbtt.plotSaveShowTimeTrace()
        else:
            print("No blobs plotted as no blobs were detected")

        # Holes
        if len(timeTraceHoles) > 0:
            if plotAll:
                loopOver = (timeTraceHolesAvg, *timeTraceHoles)
            else:
                loopOver =  (timeTraceHolesAvg, )
            for theDict in loopOver:
                pbtt.setData(theDict, "holes")
                pbtt.plotSaveShowTimeTrace()
        else:
            print("No holes plotted as no holes were detected")
#}}}

#{{{get2DData
def get2DData(ccb, varName, mode, fluct, phiCont=False):
    #{{{docstring
    """
    Driver which collects the 2D slices.

    Exists as an own entity if combined plots should be made.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    varName : str
        Name of the variable to plot.
    mode : ["perp"|"par"|"pol"]
        The mode to collect.
    fluct : bool
        Whether or not the fluctuations will be collected.
    phiCont : bool
        If True, phi contours will be overplotted.

    Returns
    -------
    blobs2DAvg : dict
        Dictionary of the averaged blob.
        Contains the keys:
            * varName    - The 2D variable.
            * varNamePPi - The 2D variable pi away from the set zInd
                           (only when mode is "par").
            * "phi"      - The 2D variable (only when phiCont is True).
            * "phiPPi"   - The 2D variable pi away from the set zInd
                           (only when mode is "par" and phiCont is True).
            * "time"     - The corresponding time.
            * "X"        - The X-mesh.
            * "Y"        - The Y-mesh.
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
        ccb.executeCollectAndCalc2D(varName, mode, fluct, phiCont)

    return blobs2DAvg, blobs2D, holes2DAvg, holes2D
#}}}

#{{{driverPlot2DData
def driverPlot2DData(ccb            ,\
                     varName        ,\
                     mode           ,\
                     fluct          ,\
                     plotSuperKwargs,\
                     plotAll        ,\
                     phiCont=False):
    #{{{docstring
    """
    Driver which plots the 2D data.

    NOTE: mode == "par" is slower then the rest, as this requires more
          opening of files.

    Parameters
    ----------
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    varName : str
        Name of the variable to plot.
    mode : ["perp"|"par"|"pol"]
        The mode to collect.
    fluct : bool
        Whether or not the fluctuations will be collected.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    plotAll : bool
       If True: The individual blobs will be plotted.
    phiCont : bool
        If True, phi contours will be overplotted.
    """
    #}}}

    # We would like all the frames to have the same max/min
    varyMaxMin = False

    blobs2DAvg, blobs2D, holes2DAvg, holes2D =\
        get2DData(ccb, varName , mode, fluct, phiCont)

    if plotAll:
        blobsAndHoles = ((blobs2DAvg, *blobs2D), (holes2DAvg, *holes2D))
    else:
        blobsAndHoles = ((blobs2DAvg, ), (holes2DAvg, ))
    for curTuple, blobOrHole in zip(blobsAndHoles,("blobs", "holes")):
        if len(curTuple[0]) == 0:
            continue
        curPSK = plotSuperKwargs.copy()
        curPSK["blobOrHole"] = blobOrHole
        for nr, blob in enumerate(curTuple):
            if nr == 0:
                curPSK["averagedBlobOrHole"] = True
            else:
                curPSK["averagedBlobOrHole"] = False
            args = (blob      ,\
                    varName   ,\
                    varyMaxMin,\
                    fluct     ,\
                    ccb       ,\
                    phiCont   ,\
                    curPSK)
            if mode == "perp":
                plotBlob2DPerp(*args)
            elif mode == "par":
                plotBlob2DPar(*args)
            elif mode == "pol":
                plotBlob2DPol(*args)
#}}}

#{{{plotBlob2DPerp
def plotBlob2DPerp(blob      ,\
                   varName   ,\
                   varyMaxMin,\
                   fluct     ,\
                   ccb       ,\
                   phiCont   ,\
                   plotSuperKwargs):
    #{{{docstring
    """
    Performs the plotting of the blob in a 2D perp plane.

    Parameters
    -----------
    blob : dict
        Dictionary with the variables to be plotted.
        Contains the keys:
            * varName - The 2D variable.
            * "phi"   - The 2D variable (only when phiCont is True).
            * "time"  - The corresponding time.
            * "X"     - The X-mesh.
            * "Y"     - The Y-mesh.
            * pos     - The position of the fixed index
    varName : str
        Name of the variable.
    varyMaxMin : bool
        Whether or not to vary the max and min for each plot in the
        colorbar.
    fluct : bool
        Wheter or not only the fluctuations are given as an input
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    phiCont : bool
        If True, phi contours will be overplotted.
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

    if phiCont:
        import pdb; pdb.set_trace()
        tupleOfArrays = (blob["phi"],)
        phiVmax, phiVmin, phiLevels =\
                getVmaxVminLevels(plotSuperKwargs,\
                                  tupleOfArrays  ,\
                                  fluct          ,\
                                  varyMaxMin     ,\
                                  nLevels = 10   ,\
                                  )

    # Burst the sequence into single pdfs
    indices = getBurstIndices(len(blob["time"]))

    for ind in indices:
        # Re-create the plotting object inside the loop to avoid
        # problems with pickling
        p2DPerp = PlotAnim2DPerp(ccb.uc          ,\
                                 fluct    = fluct,\
                                 **plotSuperKwargs)
        p2DPerp.setContourfArguments(vmax, vmin, levels)

        p2DPerp.setPerpData(blob["X"]                   ,\
                            blob["Y"]                   ,\
                            blob[varName][ind:ind+1,...],\
                            blob["time"][ind:ind+1,...] ,\
                            blob["zPos"]                ,\
                            varName)

        if phiCont:
            p2DPerp.setContourArguments(phiVmax, phiVmin, phiLevels)
            p2DPerp.setPhiData(blob["phi"][ind:ind+1,...])

        p2DPerp.plotAndSavePerpPlane()

    # Make the animation
    p2DPerp.setPerpData(blob["X"],\
                        blob["Y"],\
                        blob[varName],\
                        blob["time"],\
                        blob["zPos"],\
                        varName)
    if phiCont:
        p2DPerp.setPhiData(blob["phi"])

    p2DPerp.plotAndSavePerpPlane()
#}}}

#{{{plotBlob2DPar
def plotBlob2DPar(blob      ,\
                  varName   ,\
                  varyMaxMin,\
                  fluct     ,\
                  ccb       ,\
                  phiCont   ,\
                  plotSuperKwargs):
    #{{{docstring
    """
    Performs the plotting of the blob in a 2D par plane.

    Parameters
    -----------
    blob : dict
        Dictionary with the variables to be plotted.
        Contains the keys:
            * varName    - The 2D variable.
            * varNamePPi - The 2D variable pi away from the set zInd
            * "phi"      - The 2D variable (only when phiCont is True).
            * "phiPPi"   - The 2D variable pi away from the set zInd
                           (only when phiCont is True).
            * "time"     - The corresponding time.
            * "X"        - The X-mesh.
            * "Y"        - The Y-mesh.
            * pos        - The position of the fixed index
    varName : str
        Name of the variable.
    varyMaxMin : bool
        Whether or not to vary the max and min for each plot in the
        colorbar.
    fluct : bool
        Wheter or not only the fluctuations are given as an input
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    phiCont : bool
        If True, phi contours will be overplotted.
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

    if phiCont:
        tupleOfArrays = (blob["phi"],blob["phiPPi"])
        phiVmax, phiVmin, phiLevels =\
                getVmaxVminLevels(plotSuperKwargs,\
                                  tupleOfArrays  ,\
                                  fluct          ,\
                                  varyMaxMin     ,\
                                  nLevels = 10   ,\
                                  )

    # Burst the sequence into single pdfs
    indices = getBurstIndices(len(blob["time"]))

    for ind in indices:
        # Re-create the plotting object inside the loop to avoid
        # problems with pickling
        p2DPar = PlotAnim2DPar(ccb.uc          ,\
                               fluct    = fluct,\
                               **plotSuperKwargs)
        p2DPar.setContourfArguments(vmax, vmin, levels)


        p2DPar.setParData(blob["X"]                         ,\
                          blob["Y"]                         ,\
                          blob[varName][ind:ind+1,...]      ,\
                          blob[varName+"PPi"][ind:ind+1,...],\
                          blob["time"][ind:ind+1,...]       ,\
                          blob["thetaPos"]                  ,\
                          varName)
        if phiCont:
            p2DPol.setContourArguments(phiVmax, phiVmin, phiLevels)
            p2DPol.setPhiData   (blob["phi"][ind:ind+1,...])
            p2DPol.setPhiParData(blob["phiPPi"][ind:ind+1,...])

        p2DPar.plotAndSaveParPlane()

    # Make the animation
    p2DPar.setParData(blob["X"]          ,\
                      blob["Y"]          ,\
                      blob[varName]      ,\
                      blob[varName+"PPi"],\
                      blob["time"]       ,\
                      blob["thetaPos"]   ,\
                      varName)
    if phiCont:
        p2DPerp.setPhiData   (blob["phi"])
        p2DPerp.setPhiParData(blob["phiPPi"])

    p2DPar.plotAndSaveParPlane()
#}}}

#{{{plotBlob2DPol
def plotBlob2DPol(blob          ,\
                  varName       ,\
                  varyMaxMin    ,\
                  fluct         ,\
                  ccb           ,\
                  phiCont       ,\
                  plotSuperKwargs):
    #{{{docstring
    """
    Performs the plotting of the blob in a 2D pol plane.

    Polameters
    -----------
    blob : dict
        Dictionary with the variables to be plotted.
        Contains the keys:
            * varName - The 2D variable.
            * "phi"   - The 2D variable (only when phiCont is True).
            * "time"  - The corresponding time.
            * "X"     - The X-mesh.
            * "Y"     - The Y-mesh.
            * pos     - The position of the fixed index
    varName : str
        Name of the variable.
    varyMaxMin : bool
        Whether or not to vary the max and min for each plot in the
        colorbar.
    fluct : bool
        Wheter or not only the fluctuations are given as an input
    cbb : CollectAndCalcBlobs
        The initialized CollectAndCalcBlobs object.
    phiCont : bool
        If True, phi contours will be overplotted.
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

    if phiCont:
        tupleOfArrays = (blob["phi"],blob["phiPPi"])
        phiVmax, phiVmin, phiLevels =\
                getVmaxVminLevels(plotSuperKwargs,\
                                  tupleOfArrays  ,\
                                  fluct          ,\
                                  varyMaxMin     ,\
                                  nLevels = 10   ,\
                                  )

    # Burst the sequence into single pdfs
    indices = getBurstIndices(len(blob["time"]))

    for ind in indices:
        # Re-create the plotting object inside the loop to avoid
        # problems with pickling
        p2DPol = PlotAnim2DPol(ccb.uc          ,\
                               fluct    = fluct,\
                               **plotSuperKwargs)
        p2DPol.setContourfArguments(vmax, vmin, levels)

        p2DPol.setPolData(blob["X"]                   ,\
                          blob["Y"]                   ,\
                          blob[varName][ind:ind+1,...],\
                          blob["time"][ind:ind+1,...] ,\
                          blob["rhoPos"]              ,\
                          varName)
        if phiCont:
            p2DPol.setContourArguments(phiVmax, phiVmin, phiLevels)
            p2DPol.setPhiData(blob["phi"][ind:ind+1,...])

        p2DPol.plotAndSavePolPlane()

    # Make the animation
    p2DPol.setPolData(blob["X"],\
                      blob["Y"],\
                      blob[varName],\
                      blob["time"],\
                      blob["rhoPos"],\
                      varName)
    if phiCont:
        p2DPerp.setPhiData(blob["phi"])
    p2DPol.plotAndSavePolPlane()
#}}}

# FIXME: Could add padding pct to focus the plots closer to tau=0
#{{{getBurstIndices
def getBurstIndices(theLen, frames=9):
    #{{{docstring
    """
    Returns equally spaced indices.

    Paramters
    ---------
    theLen : int
        Number of indices in the array
    frames : int
        Number of output indices

    Returns
    -------
    indices : tuple
        Tuple of equally spaced indices.
    """
    #}}}

    # Minus 1 as the indices start from 0
    endInd = theLen-1
    # Minus 1 as we have one less segment than number
    remainder = endInd%(frames-1)
    # SUbtract the last point if the remainder is not even
    even = (remainder%2 == 0)
    if not(even):
        endInd -= 1
        remainder = endInd%(frames-1)
    # Set the padding
    pad = remainder/2

    indices = tuple(int(i) for i in np.linspace(pad, endInd-pad, frames))

    return indices
#}}}

#{{{DriverBlobs
class DriverBlobs(DriverSuperClass):
    """
    Class for driving of the plotting of the blobs.
    """
    #{{{Constructor
    def __init__(self             ,\
                 dmp_folders      ,\
                 slices           ,\
                 pctPadding       ,\
                 convertToPhysical,\
                 plotSuperKwargs  ,\
                 varName   = "n"  ,\
                 condition = 3    ,\
                 phiCont   = False,\
                 normed    = False,\
                 plotAll   = False,\
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
        varName : str
            Variable to check.
            NOTE: Condition is still on the radial density flux.
        phiCont : bool
            If True, phi contours will be overplotted.
        condition : float
            The condition in the conditional average will be set to
            flux.std()*condition
        normed : bool
            Wheter or not to norm the histogram
        plotAll : bool
           If True: All the individual frames making up the average will be
           plotted in the 2D plot.
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
        self._varName           = varName
        self._condition         = condition
        self._normed            = normed
        self._phiCont           = phiCont
        self._plotAll           = plotAll

        # Set the plot type
        plotType = os.path.join("blobs",str(condition))

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :plotType})
        self._plotSuperKwargs = plotSuperKwargs

        # Make a plotter object in order to get the picklePath
        uc = UnitsConverter(dmp_folders[0])
        tmp = PlotSuperClass(uc, **plotSuperKwargs)
        picklePath = tmp.getSavePath()

        # Prepare the blobs
        self._ccb =\
            prepareBlobs(self._collectPaths     ,\
                         slices                 ,\
                         pctPadding             ,\
                         convertToPhysical      ,\
                         condition = condition  ,\
                         picklePath = picklePath,\
                         )
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
        self.driverRadialFlux()
        modes = ("perp", "par", "pol")
        for mode in modes:
            self.setMode(mode)
            for b in (True, False):
                self.setFluct(b)
                self.driverPlot2DData()
    #}}}

    #{{{driverRadialFlux
    def driverRadialFlux(self):
        #{{{docstring
        """
        Wrapper to driverRadialFlux
        """
        #}}}
        args = (\
                self._ccb            ,\
                self._plotSuperKwargs,\
               )
        if self._useMultiProcess:
            processes = Process(target = driverRadialFlux,\
                                args   = args             )
            processes.start()
        else:
            driverRadialFlux(*args)
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
        if self._useMultiProcess:
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

        args  = (self._ccb, self._plotSuperKwargs, self._plotAll)
        if self._useMultiProcess:
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
                 self._varName        ,\
                 self._mode           ,\
                 self._fluct          ,\
                 self._plotSuperKwargs,\
                 self._plotAll        ,\
                 self._phiCont        ,\
                )
        if self._useMultiProcess:
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
        implemented = ("perp" ,"par", "pol")
        if not(mode in implemented):
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
