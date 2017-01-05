#!/usr/bin/env python

"""
Contains drivers for the radial flux
"""

from ..superClasses import DriverPointsSuperClass
from ..timeTrace import CollectAndCalcTimeTrace
from .collectAndCalcRadialFlux import CollectAndCalcRadialFlux
from .plotRadialFlux import PlotRadialFlux
from multiprocessing import Process

#{{{driverRadialFlux
def driverRadialFlux(collectPaths     ,\
                     varName          ,\
                     convertToPhysical,\
                     mode             ,\
                     indicesArgs      ,\
                     indicesKwargs    ,\
                     plotSuperKwargs  ,\
                    ):
    #{{{docstring
    """
    Driver for plotting radial fluxes.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    varName : str
        The variable name which will be used.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    mode : ["normal"|"fluct"]
        If mode is "normal" the raw data is given as an output.
        If mode is "fluct" the fluctuations are given as an output.
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
    cctt = CollectAndCalcTimeTrace(collectPaths                         ,\
                                   mode              = mode             ,\
                                   convertToPhysical = convertToPhysical,\
                                  )

    # Set the slice
    cctt.setIndices(*indicesArgs, **indicesKwargs)

    # Set name
    cctt.setVarName(varName)

    # Execute the collection
    tt = cctt.executeCollectAndCalc()
    tt = cctt.convertTo1D(tt)

    # Get the radial flux
    ccrf = CollectAndCalcRadialFlux(collectPaths, cctt.getSlices(),\
                                    mode, cctt.getDh(), cctt.convertToPhysical)
    radialExBTraces = ccrf.getRadialExBTrace()
    radialFlux      = ccrf.calcRadialFlux(tt, radialExBTraces)

    # Plot
    ptt = PlotRadialFlux(cctt.uc         ,\
                        **plotSuperKwargs)
    ptt.setData(radialFlux, mode)
    ptt.plotSaveShowRadialFlux()
#}}}

#{{{DriverRadialFlux
class DriverRadialFlux(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the radial fluxes.
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
        plotSuperKwargs.update({"plotType"   :"radialFluxes"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverRadialFlux
    def driverRadialFlux(self):
        #{{{docstring
        """
        Wrapper to driverRadialFlux
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
            processes = Process(target = driverRadialFlux, args = args)
            processes.start()
        else:
            driverRadialFlux(*args)
    #}}}
#}}}
