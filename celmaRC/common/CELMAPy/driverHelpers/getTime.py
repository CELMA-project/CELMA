#!/usr/bin/env python

"""
Contains the getTime function
"""

#{{{getTime
def getTime(depth = "second"):
    #{{{docstring
    """
    Gets the current time, and returns it as a string.

    Parameters
    ----------
    depth : ["hour" | "minute" | "second"]
        String giving the temporal accuracy of the output string.
    Returns
    -------
    nowStr : str
        The string containing the current time
    """
    #}}}
    now = datetime.datetime.now()
    nowStr = "{}-{:02d}-{:02d}".format(now.year, now.month, now.day)

    if depth == "hour" or depth == "minute" or depth == "second":
        nowStr += "-{:02d}".format(now.hour)
    if depth == "minute" or depth == "second":
        nowStr += "-{:02d}".format(now.minute)
    if depth == "second":
        nowStr += "-{:02d}".format(now.second)

    PlotSuperClass._time = nowStr
    return nowStr
#}}}
