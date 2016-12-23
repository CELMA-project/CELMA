#!/usr/bin/env python

"""
Contains the super class driver for fields 1D and fields 2D
"""

from .collectAndCalcSuperClass import CollectAndCalcSuperClass
from .plotSuperClass import PlotSuperClass

#{{{DriverPlotFieldsSuperClass
class DriverPlotFieldsSuperClass(CollectAndCalcSuperClass, PlotSuperClass):
    """
    The parent driver of 1D and 2D field plotting
    """

    #{{{Constructor
    def __init__(self           ,\
                 dmp_folders    ,\
                 xguards = False,\
                 yguards = False,\
                 xSlice  = None ,\
                 ySlice  = None ,\
                 zSlice  = None ,\
                 tSlice  = None ,\
                 xInd    = None ,\
                 yInd    = None ,\
                 zInd    = None ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:
            * Calls the parent classes
            * Sets the common memberdata

        Parameters
        ----------
        dmp_folders: tuple
            This is the output dmp_folder from bout_runners.
            Typically, these are the folders in a given scan
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        xSlice : slice
            How the data will be sliced in x.
        ySlice : slice
            How the data will be sliced in y.
        zSlice : slice
            How the data will be sliced in z.
        tSlice : slice
            How the data will be sliced in t.
        xInd : int
            Fixed rho index.
        yInd : int
            Fixed z index.
        zInd : int
            Fixed theta index.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructors of the parent classes
        # Preparing collect and plot kwargs
        collectKwargs = {}
        popKeys = ("collectPaths", "convertToPhysical")
        for key in popKeys:
            collectKwargs[key] = kwargs.pop(key)
        plotKwargs = kwargs
        plotKwargs["dmp_folders"] = dmp_folders
        CollectAndCalcSuperClass.__init__(self, dmp_folders, **collectKwargs)
        PlotSuperClass.__init__(self, self.uc, **plotKwargs)

        # Set the member data
        self._xguards = xguards
        self._yguards = yguards
        self._xSlice  = xSlice
        self._ySlice  = ySlice
        self._zSlice  = zSlice
        self._tSlice  = tSlice
        self._xInd    = xInd
        self._yInd    = yInd
        self._zInd    = zInd
    #}}}
#}}}
