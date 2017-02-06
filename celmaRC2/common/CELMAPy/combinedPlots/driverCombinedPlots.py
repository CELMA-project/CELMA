#!/usr/bin/env python

"""
Contains the driver for the combined plots
"""

from ..timeTrace import getTimeTrace
from ..radialFlux import getRadialFlux
from ..PDF import getPDF
from ..PSD import get1DPSD
from ..superClasses import DriverPointsSuperClass
from .plotCombinedPlots import PlotCombinedPlots
from multiprocessing import Process

#{{{driverCombinedPlots
def driverCombinedPlots(collectPaths     ,\
                        steadyStatePath  ,\
                        varName          ,\
                        convertToPhysical,\
                        mode             ,\
                        yInd             ,\
                        zInd             ,\
                        tSlice           ,\
                        plotSuperKwargs  ,\
                       ):
    #{{{docstring
    """
    Driver for plotting time traces.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    steadyStatePath : str
        The path to the steady state.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
    yInd : int
        Parallel fixed index.
    zInd : int
        Poloidal fixed index.
    tSlice : [None|slice]
        Time slice.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    indicesArgs   = (None, yInd, zInd)
    indicesKwargs = {\
                     "tSlice"          : tSlice         ,\
                     "nPoints"         : 3              ,\
                     "equallySpace"    : "x"            ,\
                     "steadyStatePath" : steadyStatePath,\
                    }

    args = (collectPaths     ,\
            varName          ,\
            convertToPhysical,\
            mode             ,\
            indicesArgs      ,\
            indicesKwargs    ,\
           )

    tt , uc = getTimeTrace(*args)
    rf , _  = getRadialFlux(*args)
    PDF, _  = getPDF(*args)
    PSD, _  = get1DPSD(*args)

    # Plot
    ptt = PlotCombinedPlots(uc, **plotSuperKwargs)
    ptt.setData(tt, rf, PDF, PSD, mode)
    ptt.plotSaveShowCombinedPlots()
#}}}

#{{{DriverCombinedPlots
class DriverCombinedPlots(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the combined plots.
    """

    #{{{Constructor
    def __init__(self                       ,\
                 dmp_folders                ,\
                 steadyStatePath            ,\
                 yInd                       ,\
                 zInd                       ,\
                 tSlice                     ,\
                 plotSuperKwargs            ,\
                 varName           = "n"    ,\
                 mode              = "fluct",\
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
        steadyStatePath : str
            The path to the steady state.
        yInd : int
            Parallel fixed index.
        zInd : int
            Poloidal fixed index.
        tSlice : [None|slice]
            Time slice.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        varName : str
            Name of variable to collect and plot
        mode : ["normal"|"fluct"]
            If mode is "normal" the raw data is given as an output.
            If mode is "fluct" the fluctuations are given as an output.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._varName         = varName
        self._mode            = mode
        self._steadyStatePath = steadyStatePath
        self._yInd            = yInd
        self._zInd            = zInd
        self._tSlice          = tSlice

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"combinedPlots"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverCombinedPlots
    def driverCombinedPlots(self):
        #{{{docstring
        """
        Wrapper to driverCombinedPlots
        """
        #}}}
        args = (\
                self._collectPaths    ,\
                self._steadyStatePath ,\
                self._varName         ,\
                self.convertToPhysical,\
                self._mode            ,\
                self._yInd            ,\
                self._zInd            ,\
                self._tSlice          ,\
                self._plotSuperKwargs ,\
               )
        if self._useMultiProcess:
            processes = Process(target = driverCombinedPlots, args = args)
            processes.start()
        else:
            driverCombinedPlots(*args)
    #}}}
#}}}
