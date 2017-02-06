#!/usr/bin/env python

"""
Contains single driver and driver class for the performance
"""

from ..superClasses import DriverSuperClass
from .collectAndCalcPerformance import CollectAndCalcPerformance
from .plotPerformance import PlotPerformance
from multiprocessing import Process

#{{{driverPerformance
def driverPerformance(collectPaths      ,\
                      convertToPhysical ,\
                      plotSuperKwargs   ,\
                      allFolders = False,\
                     ):
    #{{{docstring
    """
    Driver for plotting performance.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    allFolders : bool
        If true, the fileName will be tagged with allFolders.
    """
    #}}}

    # Create collect object
    ccp = CollectAndCalcPerformance(collectPaths                         ,\
                                    convertToPhysical = convertToPhysical,\
                                   )

    # Execute the collection
    perform = ccp.executeCollectAndCalc()

    ptt = PlotPerformance(ccp.uc, **plotSuperKwargs)
    ptt.setData(perform, allFolders)
    ptt.plotSaveShowPerformance()
#}}}

#{{{DriverPerformance
class DriverPerformance(DriverSuperClass):
    """
    Class for driving of the plotting of the performance.
    """

    #{{{Constructor
    def __init__(self              ,\
                 dmp_folders       ,\
                 convertToPhysical ,\
                 plotSuperKwargs   ,\
                 allFolders = False,\
                 **kwargs          ):
        #{{{docstring
        """
        This constructor:
            * Calls the parent class
            * Sets the member data
            * Updates the plotSuperKwargs

        Parameters
        ----------
        dmp_folders : tuple
            Tuple of the dmp_folder (output from bout_runners).
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        allFolders : bool
            If true, the fileName will be tagged with allFolders.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self.convertToPhysical = convertToPhysical
        self._allFolders       = allFolders

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"performance"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverPerformance
    def driverPerformance(self):
        #{{{docstring
        """
        Wrapper to driverPerformance
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self.convertToPhysical,\
                 self._plotSuperKwargs ,\
                )
        kwargs = {"allFolders":self._allFolders}
        if self._useMultiProcess:
            processes =\
                    Process(target = driverPerformance,\
                            args = args, kwargs=kwargs)
            processes.start()
        else:
            driverPerformance(*args, **kwargs)
    #}}}
#}}}
