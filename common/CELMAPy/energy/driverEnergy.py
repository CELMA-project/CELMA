#!/usr/bin/env python

"""
Contains single driver and driver class for the energies
"""

from ..superClasses import DriverSuperClass
from .collectAndCalcEnergy import CollectAndCalcEnergy
from .plotEnergy import PlotEnergy
from multiprocessing import Process

#{{{driverEnergy
def driverEnergy(collectPaths     ,\
                 convertToPhysical,\
                 tSlice           ,\
                 plotSuperKwargs  ,\
                ):
    #{{{docstring
    """
    Driver for plotting energies.

    Parameters
    ----------
    collectPaths : tuple
        Paths to collect from.
        The corresponind 't_array' of the paths must be in ascending order.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    tSlice : slice
        How to slice the time.
    plotSuperKwargs : dict
        Keyword arguments for the plot super class.
    """
    #}}}

    # Create collect object
    cce = CollectAndCalcEnergy(collectPaths                         ,\
                               convertToPhysical = convertToPhysical,\
                               )

    # Set the slice
    cce.setTSlice(tSlice)

    # Execute the collection
    e = cce.executeCollectAndCalc()

    # Plot
    pe = PlotEnergy(cce.uc          ,\
                    **plotSuperKwargs)
    pe.setData(e, timeAx = True)
    pe.plotSaveShowEnergy()
    pe.setData(e, timeAx = False)
    pe.plotSaveShowEnergy()
#}}}

#{{{DriverEnergy
class DriverEnergy(DriverSuperClass):
    """
    Class for driving of the plotting of the time traces.
    """

    #{{{Constructor
    def __init__(self                    ,\
                 dmp_folders             ,\
                 tSlice                  ,\
                 plotSuperKwargs         ,\
                 convertToPhysical = True,\
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
        plotSuperKwargs : dict
            Keyword arguments for the plot super class.
        tSlice : slice
            How to slice the time.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        **kwargs : keyword arguments
            See parent class for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(dmp_folders, **kwargs)

        # Set the member data
        self.convertToPhysical = convertToPhysical
        self._tSlice = tSlice

        # Update the plotSuperKwargs dict
        plotSuperKwargs.update({"dmp_folders":dmp_folders})
        plotSuperKwargs.update({"plotType"   :"energies"})
        self._plotSuperKwargs = plotSuperKwargs
    #}}}

    #{{{driverEnergy
    def driverEnergy(self):
        #{{{docstring
        """
        Wrapper to driverEnergy
        """
        #}}}
        args =  (\
                 self._collectPaths    ,\
                 self.convertToPhysical,\
                 self._tSlice          ,\
                 self._plotSuperKwargs ,\
                )
        if self._useMultiProcess:
            processes = Process(target = driverEnergy, args = args)
            processes.start()
        else:
            driverEnergy(*args)
    #}}}
#}}}
