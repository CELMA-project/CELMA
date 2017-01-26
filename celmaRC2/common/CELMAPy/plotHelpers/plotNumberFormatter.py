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
    # Special case if 0.000x or 0.00x
    if "0.000" in tickString:
        tickString = tickString.replace("0.000", "")
        if tickString[0] == "-":
            tickString = "{}.{}e-03".format(tickString[0:2], tickString[2:])
        else:
            tickString = "{}.{}e-03".format(tickString[0:1], tickString[1:])
    if "0.00" in tickString:
        tickString = tickString.replace("0.00", "")
        if tickString[0] == "-":
            tickString = "{}.{}e-02".format(tickString[0:2], tickString[2:])
        else:
            tickString = "{}.{}e-02".format(tickString[0:1], tickString[1:])
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
