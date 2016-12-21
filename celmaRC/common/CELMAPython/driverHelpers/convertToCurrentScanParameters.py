#!/usr/bin/env python

"""
Contains the convertToCurrentScanParameters function
"""

import re

#{{{convertToCurrentScanParameters
def convertToCurrentScanParameters(dmpFolder, aScanPath, scanParameters):
    #{{{docstring
    """
    Function which converts a path belonging to one paths in a scan
    to the path belonging to the current scan.

    The function obtains the current scan parameters from dmpFolder, and
    inserts the current scan parameters into aScanPath (the function
    input which is one of the paths belonging to the scan).

    Parameters
    ----------
    dmpFolder : str
        The dump folder under consideration.
    aScanPath : str
        One of the paths from the simulations.
    scanParameters : tuple
        A tuple of the parameters which are used for the scan (these are
        present in the path, i.e. aScanPath)

    Returns
    -------
    scanPath : str
        aScanPath converted to the scan parameters of the current run.
    """
    #}}}

    # Make a template string of aScanPath
    scanPathTemplate = aScanPath
    for scanParameter in scanParameters:
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

    # Get the values from the current dmpFolder
    values = {}
    for scanParameter in scanParameters:
        hits = [m.start() for m in \
                re.finditer(scanParameter, dmpFolder)]
        # Choose the first hit to get the value from (again we assume
        # that the value does not contain a _)
        value_start = hits[0] + len(scanParameter) + 1
        # Here we assume that the value is not separated by an
        # underscore
        values[scanParameter] = dmpFolder[value_start:].split("_")[0]

    # Insert the values
    scanPath = scanPathTemplate.format(values)

    return scanPath
#}}}
