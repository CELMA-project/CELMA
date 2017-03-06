#!/usr/bin/env python

"""
Contains single driver and driver class for the fourier modes
"""

from ..superClasses import DriverPointsSuperClass
from .collectAndCalcFourierModes import CollectAndCalcFourierModes
from .plotFourierModes import PlotFourierModes
from multiprocessing import Process

#{{{driverFourierModes
def driverFourierModes(collectPaths     ,\
                       varName          ,\
                       convertToPhysical,\
                       nModes           ,\
                       indicesArgs      ,\
                       indicesKwargs    ,\
                       plotSuperKwargs  ,\
                       ):
    #{{{docstring
    """
    Driver for plotting fourier modes.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    nModes : int
        Number of modes to plot.
    indicesArgs : tuple
        Contains xInd, yInd and zInd.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    indicesKwargs : dict
        Contains tslice, nPoints, equallySpace and steadyStatePath.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    # Create collect object
    ccfm = CollectAndCalcFourierModes(collectPaths                         ,\
                                      convertToPhysical = convertToPhysical,\
                                     )

    # Set the slice
    ccfm.setIndices(*indicesArgs, **indicesKwargs)

    # Set name
    ccfm.setVarName(varName)

    # Execute the collection
    fm = ccfm.executeCollectAndCalc()
    fm = ccfm.convertTo2D(fm)
    fm = ccfm.calcMagnitude(fm)

    # Plot
    pfm = PlotFourierModes(ccfm.uc         ,\
                           **plotSuperKwargs)
    pfm.setData(fm, nModes, timeAx = True)
    pfm.plotSaveShowFourierModes()
    pfm.setData(fm, nModes, timeAx = False)
    pfm.plotSaveShowFourierModes()
#}}}

#{{{DriverFourierModes
class DriverFourierModes(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the fourier modes.
    """

    #{{{Constructor
    def __init__(self                   ,\
                 dmp_folders            ,\
                 indicesArgs            ,\
                 indicesKwargs          ,\
                 plotSuperKwargs        ,\
                 varName           = "n",\
                 nModes            = 7  ,\
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
        indicesArgs : tuple
            Contains xInd, yInd and zInd.
            See CollectAndCalcPointsSuperClass.setIndices for details.
        indicesKwargs : dict
            Contains tslice, nPoints, equallySpace and steadyStatePath.
            See CollectAndCalcPointsSuperClass.setIndices for details.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        varName : str
            Name of variable to collect and plot
        nModes : int
            Number of modes to plot.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._varName       = varName
        self._nModes        = nModes
        self._indicesArgs   = indicesArgs
        self._indicesKwargs = indicesKwargs

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"fourierModes"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverFourierMode
    def driverFourierMode(self):
        #{{{docstring
        """
        Wrapper to driverFourierMode
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._nModes          ,\
                 self._indicesArgs     ,\
                 self._indicesKwargs   ,\
                 self._plotSuperKwargs ,\
                )
        if self._useMultiProcess:
            processes = Process(target = driverFourierModes, args = args)
            processes.start()
        else:
            driverFourierModes(*args)
    #}}}
#}}}
