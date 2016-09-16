#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .postProcessorDriver import PostProcessorDriver

#{{{StatsAndSignalsDrivers
class StatsAndSignalsDrivers(PostProcessorDriver):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self          ,\
                 *args         ,\
                 paths   = None,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets additional member data.

        Parameters
        ----------
        *args : positional arguments
            See the constructor of PostProcessorDriver for details.
        paths : list
            What folders to make a collective collect from
        **kwargs : keyword arguments
            See the constructor of PostProcessorDriver for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._paths = paths
        self._pltSize = (10, 7)
    #}}}
#}}}
