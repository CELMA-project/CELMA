#!/usr/bin/env python

"""
Contains drivers for the probes
"""

import re
from .postProcessorDriver import PostProcessorDriver

#{{{StatsAndSignalsDrivers
class StatsAndSignalsDrivers(PostProcessorDriver):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self                   ,\
                 *args                  ,\
                 paths           = None ,\
                 scanParameters  = False,\
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
        scanParameters : [None|list]
            List of parameters changed in the scan. If this is not None,
            the paths will be converted to the current scan parameters
            by calling calling convertToCurrentScanParameters.
        **kwargs : keyword arguments
            See the constructor of PostProcessorDriver for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._paths          = paths
        self._scanParameters = scanParameters
        self._pltSize        = (12, 9)

        # Convert the paths
        if self._scanParameters is not None:
            self._paths = [self._convertToCurrentScanParameters(path)
                           for path in paths]
    #}}}

    #{{{_convertToCurrentScanParameters
    def _convertToCurrentScanParameters(self, aScanPath):
        """
        Function which converts a path belonging to one paths in a scan
        to the path belonging to the current scan.

        The function obtains the current scan parameters from self._path
        (the dmp_folder given from bout_runners), and inserts the
        current scan parameters into aScanPath (the function input which is
        one of the paths belonging to the scan).

        Parameters
        ----------
        aScanPath : str
            One of the paths from the simulations.

        Returns
        -------
        scanPath : str
            aScanPath converted to the scan parameters of the current run.
        """
        # Make a template string of aScanPath
        scanPathTemplate = aScanPath
        for scanParameter in self._scanParameters:
            hits = [m.start() for m in \
                    re.finditer(scanParameter, scanPathTemplate)]
            while(len(hits) > 0):
                # Replace the values with {}
                # The value is separated from the value by 1 character
                value_start = hits[0] + len(scanParameter) + 1
                # Here we assume that the value is not separated by an
                # underscore
                value_len = len(scanPathTemplate[value_start:].split("_")[0])
                value_end = value_start + value_len
                # Replace the values with {}
                scanPathTemplate =\
                    "{}{{0[{}]}}{}".format(\
                        scanPathTemplate[:value_start],\
                        scanParameter,\
                        scanPathTemplate[value_end:])
                # Update hits
                hits.remove(hits[0])

        # Get the values from the current self._path
        values = {}
        for scanParameter in self._scanParameters:
            hits = [m.start() for m in \
                    re.finditer(scanParameter, self._path)]
            # Choose the first hit to get the value from (again we assume
            # that the value does not contain a _)
            value_start = hits[0] + len(scanParameter) + 1
            # Here we assume that the value is not separated by an
            # underscore
            values[scanParameter] = self._path[value_start:].split("_")[0]

        # Insert the values
        scanPath = scanPathTemplate.format(values)

        return scanPath
    #}}}
#}}}
