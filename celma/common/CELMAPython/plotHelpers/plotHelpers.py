#!/usr/bin/env python

""" Collection of routines which helps the plotting """

import scipy.constants as cst

#{{{plotNumberFormatter
def plotNumberFormatter(val, pos):
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

    tickString = "${:.3g}".format(val)
    if "e+" in tickString:
        tickString = tickString.replace("e+" , r"\cdot 10^{")
        tickString = tickString.replace("e+0", r"\cdot 10^{")
        tickString += "}$"
    elif "e-" in tickString:
        tickString = tickString.replace("e-" , r"\cdot 10^{-")
        tickString = tickString.replace("e-0", r"\cdot 10^{-")
        tickString += "}$"
    else:
        tickString += "$"

    return tickString
#}}}

#{{{physicalUnitsConverter
def physicalUnitsConverter(var, varName, convertToPhysical, convDict):
    #{{{docstring
    """
    Calculates physical parameters from the normalized if
    convertToPhysical is set. Returns the units.

    Parameters
    ----------
    var : array
        The variable.
    varName : str
        Name of the variable.
    convertToPhysical : bool
        Whether or not to convert to physical units.
    convDict : dict
        Dictionary containing the values of:
            * omCI
            * rhoS
            * n0
            * Te0

    Returns
    -------
    var : array
        The variable after eventual processing.
    units : str
        The units which will be plotted (does not contain $ for
        LaTeX.
    """
    #}}}

    if convertToPhysical:
        # Calculate back to physical units
        if varName == "n":
            var *= convDict["n0"]
            units = r"\mathrm{m}^{-3}"
        elif varName == "vort":
            var *= convDict["omCI"]
            units = r"\mathrm{s}^{-1}"
        elif varName == "vortD":
            var *= convDict["omCI"]*convDict["n0"]
            units = r"\mathrm{m}^{-3}\mathrm{s}^{-1}"
        elif varName == "phi":
            var *= convDict["Te0"]/cst.e
            units = r"\mathrm{J}\mathrm{C}^{-1}"
        elif varName == "jPar":
            var *= cst.m_p*convDict["n0"]*convDict["rhoS"]*convDict["n0"]
            units = r"\mathrm{C}\mathrm{s}^{-1}"
        elif varName == "momDensPar":
            # momDensPar is divided by m_i, so we need to multiply
            # by m_i again here
            var *= cst.m_p*convDict["rhoS"]*convDict["omCI"]*convDict["n0"]
            units = r"\mathrm{kg }\mathrm{m}^{-2}\mathrm{s}^{-1}"
        elif varName == "uIPar":
            var *= convDict["rhoS"]*convDict["omCI"]
            units = r"\mathrm{m}\mathrm{s}^{-1}"
        elif varName == "uEPar":
            var *= convDict["rhoS"]*convDict["omCI"]
            units = r"\mathrm{m}\mathrm{s}^{-1}"
        elif varName == "S":
            var *= convDict["omCI"]*convDict["n0"]
            units = r"\mathrm{m}^{-3}\mathrm{s}^{-1}"
        elif varName == "t":
            var /= convDict["omCI"]
            units = r"\mathrm{s}"
        elif varName == "rho":
            var *= convDict["rhoS"]
            units = r"\mathrm{m}"
        elif varName == "z":
            var *= convDict["rhoS"]
            units = r"\mathrm{m}"
        else:
            units = " "
    else:
        # Calculate back to physical units
        if varName == "n":
            units = r"/n_0"
        elif varName == "vort":
            units = r"/\omega_{{ci}}"
        elif varName == "vortD":
            units = r"/\omega_{{ci}}n_0"
        elif varName == "phi":
            units = r" q/T_{{e,0}}"
        elif varName == "jPar":
            units = r"/n_0c_sq"
        elif varName == "momDensPar":
            # by m_i again here
            units = r"/m_in_0c_s"
        elif varName == "uIPar":
            units = r"/c_s"
        elif varName == "uEPar":
            units = r"/c_s"
        elif varName == "S":
            units = r"/\omega_{{ci}}n_0"
        elif varName == "t":
            units = r"\omega_{{ci}}"
        elif varName == "rho":
            units = r"/\rho_s"
        elif varName == "z":
            units = r"/\rho_s"
        else:
            units = " "

    return var, units
#}}}
