#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import collectEnergy, PlotEnergy

#{{{DriversEnergy
class DriversEnergy(StatsAndSignalsDrivers):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self                      ,\
                 *args                     ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See the constructor of StatsAndSignalsDrivers for details.
        **kwargs : keyword arguments
            See the constructor of StatsAndSignalsDrivers for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Placeholder for the energy
        self._energy = None
    #}}}

    #{{{collectEnergies
    def collectEnergies(self):
        """ Collect the energy from the files """
        # Create the probes
        self._energy = collectEnergy(self._paths)
    #}}}

    #{{{plotEnergy
    def plotEnergy(self):
        """ Plot the energy """

        # Calculate the probes if not already done
        if self._energy == None:
            self.collectEnergies()

        # Create the energyPlotter
        energyPlotter = PlotEnergy(\
                self._paths                                ,\
                self._energy                               ,\
                convertToPhysical = self.convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        energyPlotter.plotKinEnergy("ions")
        energyPlotter.plotKinEnergy("electrons")
        energyPlotter.plotPotEnergy()
    #}}}
#}}}
