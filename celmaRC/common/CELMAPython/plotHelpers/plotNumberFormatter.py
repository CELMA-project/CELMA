#!/usr/bin/env python

""" Contains function which make the ticks nice """

#{{{plotNumberFormatter
def plotNumberFormatter(val, pos, precision=3):
    #{{{docstring
    """
    Formatting numbers in the plot

    Parameters
    ----------
    val : float
        The value.
    pos : [None | float]
        The position (needed as input from FuncFormatter).
    """
    #}}}

    tickString = "${{:.{}g}}".format(precision).format(val)
    if "e+" in tickString:
        tickString = tickString.replace("e+0", r"\cdot 10^{")
        tickString = tickString.replace("e+" , r"\cdot 10^{")
        tickString += "}$"
    elif "e-" in tickString:
        tickString = tickString.replace("e-0", r"\cdot 10^{-")
        tickString = tickString.replace("e-" , r"\cdot 10^{-")
        tickString += "}$"
    else:
        tickString += "$"

    return tickString
#}}}
