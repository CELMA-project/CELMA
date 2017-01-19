#!/usr/bin/env python

"""
Contains drivers for the total flux
"""

from ..superClasses import DriverPointsSuperClass
from .collectAndCalcTotalFlux import CollectAndCalcTotalFlux
from .plotTotalFlux import PlotTotalFlux
from multiprocessing import Process

#{{{driverTotalFlux
def driverTotalFlux(collectPaths     ,\
                    xInd             ,\
                    yInd             ,\
                    tSlice           ,\
                    mode             ,\
                    convertToPhysical,\
                    plotSuperKwargs  ,\
                   ):
    #{{{docstring
    """
    Driver for plotting total fluxes.

    Parameters
    ----------
    collectPaths : tuple
        Tuple from where to collect
    xInd : int
        How to slice in the radial direction.
    yInd : int
        How to slice in the parallel direction.
    tSlice : [None|slice]
        How to slice in time.
    mode : ["normal"|"fluct"]
        Whether to look at fluctuations or normal data
    convertToPhysical : bool
        Whether or not to convert to physical
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    ccTF = CollectAndCalcTotalFlux(\
                                   collectPaths                         ,\
                                   xInd                                 ,\
                                   yInd                                 ,\
                                   mode              = mode             ,\
                                   convertToPhysical = convertToPhysical,\
                                  )

    intFluxes = ccTF.executeCollectAndCalc()


    # Plot
    import pdb; pdb.set_trace()
    ptf = PlotTotalFlux(ccTF.uc         ,\
                        **plotSuperKwargs)
    ptf.setData(intFluxes, mode)
    ptf.plotSaveShowTotalFlux()
#}}}

#{{{DriverTotalFlux
class DriverTotalFlux(DriverPointsSuperClass):
    """
    Class for driving of the plotting of the total fluxes.
    """

    #{{{Constructor
    def __init__(self                        ,\
                 dmp_folders                 ,\
                 xInd                        ,\
                 yInd                        ,\
                 tSlice                      ,\
                 plotSuperKwargs             ,\
                 mode              = "normal",\
                 convertToPhysical = True    ,\
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
        collectPaths : tuple
            Tuple from where to collect
        xInd : int
            How to slice in the radial direction.
        yInd : int
            How to slice in the parallel direction.
        tSlice : [None|slice]
            How to slice in time.
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        mode : ["normal"|"fluct"]
            Whether to look at fluctuations or normal data
        convertToPhysical : bool
            Whether or not to convert to physical
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self._xInd   = xInd
        self._yInd   = yInd
        self._tSlice = tSlice
        self._mode   = mode

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"totalFluxes"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverTotalFlux
    def driverTotalFlux(self):
        #{{{docstring
        """
        Wrapper to driverTotalFlux
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self._xInd            ,\
                 self._yInd            ,\
                 self._tSlice          ,\
                 self._mode            ,\
                 self.convertToPhysical,\
                 self._plotSuperKwargs ,\
                )
        if self._useSubProcess:
            processes = Process(target = driverTotalFlux, args = args)
            processes.start()
        else:
            driverTotalFlux(*args)
    #}}}
#}}}
