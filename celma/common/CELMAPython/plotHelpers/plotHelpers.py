#!/usr/bin/env python

""" Collection of routines which helps the plotting """

import scipy.constants as cst

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
        The normalization which will be plotted. Does not contain $ for
        LaTeX. An empty string is returned if convertToPhysical is True.
    units : str
        The units which will be plotted. Does not contains $ for
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
            # Generic for velocities
        elif varName == "u":
            var *= convDict["rhoS"]*convDict["omCI"]
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
        units = ""
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

#{{{getCoordinates
def getCoordinates(path                     ,\
                   xguards           = False,\
                   yguards           = False,\
                   convertToPhysical = False,\
                   ):
    #{{{docstring
    """
    Get the coordinates.

    If convertToPhysical = True, the variables in normalized units is
    returned.

    NOTE: t is collected on its own due to the possibility to slice in
          t.

    Parameters
    ----------
    path : str
        The path to collect from.
    xguards : bool
        If xguards should be included when collecting.
    yguards : bool
        If yguards should be included when collecting.
    convertToPhysical : bool
        If the physical or normalized units should be plotted.

    Returns
    -------
    The output is returned as the dictionary returnDict, which has the
    following keys:

    rho : array
        The radial positions.
    theat : array
        The poloidal positions.
    z : array
        The positions along the magnetic field.
    """
    #}}}

    #{{{rho
    dx = collect("dx"             ,\
                 path    = path   ,\
                 xguards = xguards,\
                 yguards = yguards,\
                 info    = False)
    MXG = collect("MXG"            ,\
                  path    = path   ,\
                  xguards = xguards,\
                  yguards = yguards,\
                  info    = False)

    nPoints = dx.shape[0]
    dx      = dx[0,0]

    if xguards:
        innerPoints = nPoints - 2*MXG
    else:
        innerPoints = nPoints

    # By default there is no offset in the cylinder
    # For comparision with other codes, an offset option is set
    # Read the input file
    myOpts = BOUTOptions(path)
    # Read in geom offset
    try:
        offset = eval(myOpts.geom["offset"])
        spacing = "\n"*3
        print("{0}!!!WARNING: 'offset' found in BOUT.inp, "
              "running as annulus!!!{0}".format(spacing))
        rho = offset + dx * np.array(np.arange(0.5, innerPoints))
    except KeyError:
        # This is the default
        rho = dx * np.array(np.arange(0.5, innerPoints))

    if xguards:
        # Insert the first and last grid point
        rho = np.insert(rho, 0, - 0.5*dx)
        rho = np.append(rho, rho[-1] + dx)
    #}}}

    #{{{z
    dy  = collect("dy"             ,\
                  path    = path   ,\
                  xguards = xguards,\
                  yguards = yguards,\
                  info    = False)
    MYG = collect("MYG"            ,\
                  path    = path   ,\
                  xguards = xguards,\
                  yguards = yguards,\
                  info    = False)

    nPoints  = dy.shape[1]
    dy = dy[0,0]

    if yguards:
        innerPoints = nPoints - 2*MYG
    else:
        innerPoints = nPoints

    z = dy * np.array(np.arange(0.5, innerPoints))

    if yguards:
        # Insert the first and last grid point
        z = np.insert(z, 0, - 0.5*dy)
        z = np.append(z, z[-1] + dy)
    #}}}

    #{{{theta
    dz = collect("dz"             ,\
                 path    = path   ,\
                 xguards = xguards,\
                 yguards = yguards,\
                 info    = False)
    MZ       = collect("MZ"             ,\
                       path    = path   ,\
                       xguards = xguards,\
                       yguards = yguards,\
                       info    = False)

    # Subtract the unused plane
    innerPoints = MZ - 1

    theta = dz * np.array(np.arange(0.0, innerPoints))

    # Convert to degrees
    theta * (180/np.pi)
    #}}}

    returnDict = {"rho"  :rho  ,\
                  "theta":theta,\
                  "z"    :z    ,\
                 }

    return returnDict
#}}}

#{{{getConversionDict
def getConversionDict(path, convertToPhysical):
    #{{{docstring
    """
    Get the conversion dictionary.

    Parameters
    ----------
    path : str
        The path to collect from.
    convertToPhysical : bool
        If the physical or normalized units should be plotted.

    Returns
    -------
    convDict : dictionary
        A dictionary used when converting units
    convertToPhysical : bool
        Will be set to False if the normalization units are not found.
    """
    #}}}

        convDict = {}
        if convertToPhysical:
            try:
                normalizers = ["omCI", "rhoS", "n0", "Te0"]
                for normalizer in normalizers:
                    convDict[normalizer] =\
                            collect(normalizer, path=path, info=False)

            except ValueError:
                # An OSError is thrown if the file is not found
                message = ("{0}{1}WARNING: Normalized quantities not found. "
                           "The time remains normalized".format("\n"*3,"!"*3))
                print(message)

                # Reset convertToPhysical
                convertToPhysical = False

    return convDict, convertToPhysical
#}}}

#{{{coordinatesConvertAndGetLabels
def coordinatesConvertAndGetLabels(path, convertToPhysical):
    # FIXME: YOU ARE HERE
    #{{{docstring
    """
    Get the conversion dictionary.

    Parameters
    ----------
    path : str
        The path to collect from.
    convertToPhysical : bool
        If the physical or normalized units should be plotted.

    Returns
    -------
    convDict : dictionary
        A dictionary used when converting units
    convertToPhysical : bool
        Will be set to False if the normalization units are not found.
    """
    #}}}

        convDict = {}
        if convertToPhysical:
            try:
                normalizers = ["omCI", "rhoS", "n0", "Te0"]
                for normalizer in normalizers:
                    convDict[normalizer] =\
                            collect(normalizer, path=path, info=False)

            except ValueError:
                # An OSError is thrown if the file is not found
                message = ("{0}{1}WARNING: Normalized quantities not found. "
                           "The time remains normalized".format("\n"*3,"!"*3))
                print(message)

                # Reset convertToPhysical
                convertToPhysical = False

    return convDict, convertToPhysical
#}}}
