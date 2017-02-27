#!/usr/bin/env python

"""
Contains single driver and driver class for the performance
"""

from ..superClasses import DriverSuperClass
from .collectAndCalcPerformance import CollectAndCalcPerformance
from .plotPerformance import PlotPerformance
from multiprocessing import Process

#{{{driverPerformance
def driverPerformance(collectPaths     ,\
                      convertToPhysical,\
                      mode             ,\
                      plotSuperKwargs  ,\
                      tSlice = None    ,\
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
    mode : ["init"|"expand"|"linear"|"turbulence"|"all"]
        What part of the simulation is being plotted for.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    tSlice : slice
        Use if the data should be sliced.
    """
    #}}}

    # Create collect object
    ccp = CollectAndCalcPerformance(collectPaths                         ,\
                                    convertToPhysical = convertToPhysical,\
                                   )

    # Execute the collection
    perform = ccp.executeCollectAndCalc(tSlice = tSlice)

    ptt = PlotPerformance(ccp.uc, **plotSuperKwargs)
    ptt.setData(perform, mode)
    ptt.plotSaveShowPerformance()
#}}}

#{{{DriverPerformance
class DriverPerformance(DriverSuperClass):
    """
    Class for driving of the plotting of the performance.
    """

    #{{{Constructor
    def __init__(self             ,\
                 dmp_folders      ,\
                 convertToPhysical,\
                 mode             ,\
                 plotSuperKwargs  ,\
                 tSlice = None    ,\
                 **kwargs):
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
        mode : ["init"|"expand"|"linear"|"turbulence"|"all"]
            What part of the simulation is being plotted for.
        tSlice : slice
            Use if the data should be sliced.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self.convertToPhysical = convertToPhysical
        self._mode             = mode
        self._tSlice           = tSlice

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
                 self._mode            ,\
                 self._plotSuperKwargs ,\
                )
        kwargs = {"tSlice":self._tSlice}
        if self._useMultiProcess:
            processes =\
                    Process(target = driverPerformance,\
                            args = args, kwargs=kwargs)
            processes.start()
        else:
            driverPerformance(*args, **kwargs)
    #}}}
#}}}
