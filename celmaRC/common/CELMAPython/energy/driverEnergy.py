#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import collectEnergy, PlotEnergy

#{{{DriverEnergy
class DriverEnergy(CommonPostProcessingDriver):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self,\
                 *args,\
                 paths,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets path and pltSize

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        paths  : tuple of strings
            The paths to collect from
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._paths   = paths
        self._pltSize = (12, 9)

        # Placeholder for the energy
        self._energy = None
    #}}}

    #{{{getEnergies
    def getEnergies(self):
        """Obtain the energies"""
        # Create the probes
        self._energy = calcEnergy(self._paths)
    #}}}

    #{{{plotEnergy
    def plotEnergy(self):
        """ Plot the energy """

        # Calculate the probes if not already done
        if self._energy == None:
            self.getEnergies()

        # Create the energyPlotter
        energyPlotter = PlotEnergy(\
                self._paths                                ,\
                self._energy                               ,\
                convertToPhysical = self._convertToPhysical,\
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
