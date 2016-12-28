#!/usr/bin/env python

"""
Contains helper functions used when dealing with scans
"""

#{{{getScanValue
def getScanValue(path, scanParameter):
    #{{{docstring
    """
    Returns the scan value from a path

    Parameters
    ----------
    path : [str | tuple of string]
        Path or paths referring to one particular scan value
    scanParameter : str
        String segment to search for in the path name

    Returns
    -------
    scanValue : float
        The value used in the scan.
    """
    #}}}

    # Obtain the current scan value
    if hasattr(path, "__iter__") and type(path) != str:
        # The path is a sequence
        scanValStr = path[0].split("_")
    else:
        # The path is a string
        scanValStr = path.split("_")

    # +1 as the value is immediately after the scan parameter
    scanValue = scanValStr[scanValStr.index(scanParameter)+1]

    import pdb; pdb.set_trace()
    return float(scanValue)
#}}}
