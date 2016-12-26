#!/usr/bin/env python

"""
Contains single driver and driver class for the fourier modes
"""

from ..superClasses import DriverPointsSuperClass
from .collectAndCalcFourierMode import CollectAndCalcFourierMode
from .plotFourierMode import PlotFourierMode
import os

#{{{driverFourierMode
def driverFourierMode(collectPaths     ,\
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
    ccfm = CollectAndCalcFourierMode(collectPaths                         ,\
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
    pfm = PlotFourierMode(ccfm.uc         ,\
                          **plotSuperKwargs)
    pfm.setData(fm, nModes)
    pfm.plotSaveShowFourierMode()
#}}}

#{{{DriverFourierMode
class DriverFourierMode(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the fourier modes.
    """

    #{{{Constructor
    def __init__(self                       ,\
                 dmp_folders                ,\
                 indicesArgs                ,\
                 indicesKwargs              ,\
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
        self._varName       = varName
        self._mode          = mode
        self._indicesArgs   = indicesArgs
        self._indicesKwargs = indicesKwargs

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
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
                 self._mode            ,\
                 self._indicesArgs     ,\
                 self._indicesKwargs   ,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverFourierMode, args = args)
            processes.start()
        else:
            driverFourierMode(*args)
    #}}}
#}}}
