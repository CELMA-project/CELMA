#!/usr/bin/env python

"""
Contains drivers for the blobs
"""

from ..superClasses import DriverSuperClass
from .collectAndCalcBlobs import CollectAndCalcBlobs
from .plotBlobs import PlotBlobs, PlotTemporalStats
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

    holesWaitingTime, holesPulseWidths = ccb.getWaitingTimesAndPulseWidth("holes")
    blobsWaitingTime, blobsPulseWidths = ccb.getWaitingTimesAndPulseWidth("blobs")

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

    timeTraceBlobAvg, timeTraceBlobs, timeTraceHolesAvg, timeTraceHoles =\
        ccb.executeCollectAndCalc1D()

# FIXME: Auto detect whether average or not
    ptt = PlotBlobs(ccb.uc              ,\
                    **plotSuperKwargs)
    ptt.setData(blobBinsDict, mode="foo")
    ptt.plotSaveShowBlobs()
#}}}

# FIXME:
#{{{get2DData
def get2DData(ccb, mode, fluct):
    #{{{docstring
    """
    Driver which collects the 2D slices.

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

    return blobsAvg2D, blobs2D, holesAvg2D, holes2D
#}}}

# FIXME:
#{{{driverPlot2DData
def driverPlot2DData(data2D, includeBins, plotSuperKwargs):
    #{{{docstring
    """
    Driver which plots the 2D data.

    Parameters
    ----------
    data2D : tuple
        Tuple containing blobs2DAvg, blobs2D, holes2DAvg and holes2D.
    includeBins : bool
        Whether or not the bins should be included in the plot.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}


# FIXME: Can probably use the already made 2D plotter for this
# FIXME: Just adjust the names
    ptt = PlotBlobs(ccb.uc              ,\
                    **plotSuperKwargs)
    ptt.setData(blobBinsDict, mode="foo")
    ptt.plotSaveShowBlobs()
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

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"blobs"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverBlobs
    def driverBlobs(self):
        #{{{docstring
        """
        Wrapper to driverBlobs
        """
        #}}}
        args =  (\
                 self._collectPaths     ,\
                 self._slices           ,\
                 self._pctPadding       ,\
                 self._convertToPhysical,\
                 self._plotSuperKwargs  ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverBlobs, args = args)
            processes.start()
        else:
            driverBlobs(*args)
    #}}}
#}}}
