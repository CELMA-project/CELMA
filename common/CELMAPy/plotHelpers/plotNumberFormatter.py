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
    checkForPeriod = False
    # Special case if 0.000x or 0.00x
    if "0.000" in tickString:
        checkForPeriod = True
        tickString = tickString.replace("0.000", "")
        if tickString[1] == "-":
            tickString = "{}.{}e-03".format(tickString[0:3], tickString[3:])
        else:
            tickString = "{}.{}e-03".format(tickString[0:2], tickString[2:])
    elif "0.00" in tickString:
        checkForPeriod = True
        tickString = tickString.replace("0.00", "")
        if tickString[1] == "-":
            tickString = "{}.{}e-02".format(tickString[0:3], tickString[3:])
        else:
            tickString = "{}.{}e-02".format(tickString[0:2], tickString[2:])
    if checkForPeriod:
        # The last character before e-0x should not be a period
        if len(tickString)>5 and tickString[-5] == ".":
            tickString = tickString.replace(".","")
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
