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
    normalization : str
        The normalization which will be plotted  (does not contain $ for
        LaTeX. An empty string is returned if convertToPhysical is True.
    units : str
        The units which will be plotted (does not contain $ for
        LaTeX. An empty string is returned if convertToPhysical is False.
    """
    #}}}

    if convertToPhysical:
        normalization = ""
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
        untis = ""
        # Return normalization
        if varName == "n":
            normalization = r"/n_0"
        elif varName == "vort":
            normalization = r"/\omega_{{ci}}"
        elif varName == "vortD":
            normalization = r"/\omega_{{ci}}n_0"
        elif varName == "phi":
            normalization = r" q/T_{{e,0}}"
        elif varName == "jPar":
            normalization = r"/n_0c_sq"
        elif varName == "momDensPar":
            normalization = r"/m_in_0c_s"
        elif varName == "uIPar":
            normalization = r"/c_s"
        elif varName == "uEPar":
            normalization = r"/c_s"
        elif varName == "S":
            normalization = r"/\omega_{{ci}}n_0"
        elif varName == "t":
            normalization = r"\omega_{{ci}}"
        elif varName == "rho":
            normalization = r"/\rho_s"
        elif varName == "z":
            normalization = r"/\rho_s"
        else:
            normalization = " "

    return var, normalization, units
#}}}
