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
    convertToPhysical : bool
        Whether or not to convert to physical
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    ccb = CollectAndCalcBlobs(collectPaths, slices, convertToPhysical)

    blobBinsDict, holeBinsDict, blobsAvgDict, holesAvgDict =\
        ccb.executeCollectAndCalc()

    import pdb; pdb.set_trace()

    # Plot
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
                 self._slices           ,\
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
