#!/usr/bin/env python

"""
Contains single driver and driver class for the magnitude spectrum
"""

from ..superClasses import DriverPointsSuperClass
from .collectAndCalcKThetaSpectrum import CollectAndCalcKThetaSpectrum
from .plotKThetaSpectrum import PlotKThetaSpectrum
from multiprocessing import Process

#{{{driverKThetaSpectrum
def driverKThetaSpectrum(collectPaths     ,\
                            varName          ,\
                            convertToPhysical,\
                            indicesArgs      ,\
                            indicesKwargs    ,\
                            plotSuperKwargs  ,\
                           ):
    #{{{docstring
    """
    Driver for plotting magnitude spectrum.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    indicesArgs : tuple
        Contains xInd, yInd and zInd.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    indicesKwargs : dict
        Contains tslice, nPoints, equallySpace and steadyStatePath.
        See CollectAndCalcPointsSuperClass.setIndices for details.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    includeErrorBars : bool
        Whether or not to include the error bars.
    """
    #}}}

    # Create collect object
    ccms = CollectAndCalcKThetaSpectrum(\
                                collectPaths     ,\
                                varName          ,\
                                convertToPhysical,\
                                indicesArgs      ,\
                                indicesKwargs    ,\
                                         )

    # Execute the collection
    mSpec, uc = ccms.executeCollectAndCalc()

    # Plot
    for b in (True, False):
        pfm = PlotKThetaSpectrum(uc                  ,\
                                    includeErrorBars = b,\
                                    **plotSuperKwargs)
        pfm.setData(mSpec, varName)
        pfm.plotSaveShowKThetaSpectrum()
#}}}

#{{{DriverKThetaSpectrum
class DriverKThetaSpectrum(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the magnitude spectrum.
    """

    #{{{Constructor
    def __init__(self                   ,\
                 dmp_folders            ,\
                 indicesArgs            ,\
                 indicesKwargs          ,\
                 plotSuperKwargs        ,\
                 varName           = "n",\
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
        self._indicesArgs   = indicesArgs
        self._indicesKwargs = indicesKwargs

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"kThetaSpectrum"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverKThetaSpectrum
    def driverKThetaSpectrum(self):
        #{{{docstring
        """
        Wrapper to driverFourierMode
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self._varName         ,\
                 self.convertToPhysical,\
                 self._indicesArgs     ,\
                 self._indicesKwargs   ,\
                 self._plotSuperKwargs ,\
                )
        if self._useMultiProcess:
            processes = Process(target = driverKThetaSpectrum, args = args)
            processes.start()
        else:
            driverKThetaSpectrum(*args)
    #}}}
#}}}
