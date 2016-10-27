#!/usr/bin/env python

"""
Contains drivers for plotting energy and probes
"""

from .driversEnergy import DriversEnergy
from .driversProbes import DriversProbes
from multiprocessing import Process

#{{{DriversEnergyProbes
class DriversEnergyProbes(DriversProbes, DriversEnergy):
    """
    Class which combines energy and probes plots.
    """

    #{{{Constructor
    def __init__(self             ,\
                 *args            ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See the constructor of DriversProbes and DriversEnergy for details.
        **kwargs : keyword arguments
            See the constructor of DriversProbes and DriversEnergy for details.
        """
        #}}}

        # Call the constructor of the parent class
        super(DriversEnergyProbes, self).__init__(*args, **kwargs)
    #}}}

    #{{{plotEnergyAndProbes
    def plotEnergyAndProbes(self):
        if self._useSubProcess:
            #{{{ Function call through subprocess
            # Plot the energy
            Process(\
                    target = self.plotEnergy,\
                    args   = ()             ,\
                    kwargs = {}
                   ).start()

            # Plot the probes
            Process(\
                    target = self.plotProbes,\
                    args   = ()             ,\
                    kwargs = {}
                   ).start()
            #}}}
        else:
            #{{{ Normal function call
            # Plot the energy
            self.plotEnergy()

            # Plot the probes
            self.plotProbes()
            #}}}

    #}}}
#}}}
