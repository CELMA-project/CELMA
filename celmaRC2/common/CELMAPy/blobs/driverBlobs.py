#!/usr/bin/env python

"""
Contains drivers for the blobs
"""

from ..superClasses import DriverSuperClass
from .collectAndCalcBlobs import CollectAndCalcBlobs
from .plotBlobs import PlotBlobs
from multiprocessing import Process

#{{{driverBlobs
def driverBlobs(collectPaths     ,\
                slices           ,\
                pctPadding       ,\
                convertToPhysical,\
                plotSuperKwargs  ,\
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
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    ccb = CollectAndCalcBlobs(collectPaths, slices, convertToPhysical)

    ccb.prepareCollectAndCalc()

    rf = ccb.getRadialFlux()

    holesWaitingTime, holesPulseWidths = ccb.getWaitingTimesAndPulseWidth("holes")
    blobsWaitingTime, blobsPulseWidths = ccb.getWaitingTimesAndPulseWidth("blobs")

    timeTraceBlobAvg, timeTraceBlobs, timeTraceHolesAvg, timeTraceHoles =\
        ccb.executeCollectAndCalc1D()

    perpBlobAvg, perpBlobs, perpHolesAvg, perpHoles =\
        ccb.executeCollectAndCalc2D("perp", False)

    YOU ARE HERE, TEST PAR AND POL
    import pdb; pdb.set_trace()
    a = 1

    # Plot
# FIXME: Have different kind of plots here
# First: Holes and blobs: Can use the same drivers, just add txt "blob", "hole"
# Second: bins and averages
# Third: 1D and 2D
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
