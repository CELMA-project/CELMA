#!/usr/bin/env python

"""
Contains function which calculates the energies
"""

from ..plotHelpers import PlotHelper, collectiveCollect, seqCMap, seqCMap3
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os

#{{{collectAndCalcEnergy
def collectAndCalcEnergy(paths                   ,\
                         convertToPhysical = True,\
                         tSlice            = None):
    #{{{docstring
    """
    Collects and concatenate the energies

    Parameters
    ----------
    paths : iterable of strings
        The paths to collect from. Must be in ascending order of the
        simulation time, as the variables are being concatenated
    convertToPhysical : bool
        Whether or not to convert to physical units.
    tSlice : [None|Slice}
        Whether or not to slice the time trace

    Return
    ------
    energies : dict
FIXME: Bug in monitors, needs to be fixed
FIXME: No longer support for Helmholtz like energy
        A dictionary of the energies containing the following keys
        * fluctPerpKinEE  - Fluctuation of perpendicular kinetic electron energy
        * fluctParKinEE   - Fluctuation of parallel kinetic electron energy
        * fluctSumKinEE   - Fluctuation of sum of the kinetic electron energy
        * fluctPerpKinEI  - Fluctuation of perpendicular kinetic ion energy
        * fluctParKinEI   - Fluctuation of parallel kinetic ion energy
        * fluctSumKinEI   - Fluctuation of sum of the kinetic ion energy
        * fluctAvgPotEE   - Theta avg of ppotential electron energy (from nT)
        * polAvgPerpKinEE - Theta avg perpendicular kinetic electron energy
        * polAvgParKinEE  - Theta avg parallel kinetic electron energy
        * polAvgSumKinEE  - Theta avg sum of the kinetic electron energy
        * polAvgPerpKinEI - Theta avg perpendicular kinetic ion energy
        * polAvgParKinEI  - Theta avg parallel kinetic ion energy
        * polAvgSumKinEI  - Theta avg sum of the kinetic ion energy
        * polAvgPotEE     - Theta avg of ppotential electron energy (from nT)
        * totPerpKinEE    - Total perpendicular kinetic electron energy
        * totParKinEE     - Total parallel kinetic electron energy
        * totSumKinEE     - Total sum of the kinetic electron energy
        * totPerpKinEI    - Total perpendicular kinetic ion energy
        * totParKinEI     - Total parallel kinetic ion energy
        * totSumKinEI     - Total sum of the kinetic ion energy
        * t               - Time
    """
    #}}}

# FIXME: t_array can be collected differently
"t_array"        ,\
    varStrings= (\
                 "polAvgPerpKinEE",\
                 "polAvgParKinEE" ,\
                 "polAvgSumKinEE" ,\
                 "polAvgPerpKinEI",\
                 "polAvgParKinEI" ,\
                 "polAvgSumKinEI" ,\
                 "fluctPerpKinEE" ,\
                 "fluctParKinEE"  ,\
                 "fluctSumKinEE"  ,\
                 "fluctPerpKinEI" ,\
                 "fluctParKinEI"  ,\
                 "fluctSumKinEI"  ,\
                 "polAvgPotEE"    ,\
                 "fluctAvgPotEE"  ,\
                 "t_array"        ,\
                )

        if tSlice is not None:
            tStart = tSlice[tCounter].start
            tEnd   = tSlice[tCounter].stop
            tCounter += 1
        else:
            tStart = None
            tEnd   = None
# FIXME: Waiting for a bug here, when occurs, fallback to python implementation
    energies = collectiveCollect(paths, varStrings)

    # Change keyname of t_array to t
    energies["t"] = energies.pop("t_array")

    if tSlice is not None:
        # Slice the variables with the step
        # Make a new slice as the collect dealt with the start and
        # the stop of the slice
        newSlice = slice(None, None, tSlice.step)
        for key in energies.keys():
            energies[key] = energies[key][newSlice]

# FIXME: Fix this!!!
#    firstCharToLower = lambda s: s[:1].lower() + s[1:] if s else ""
#    firstCharToUpper = lambda s: s[:1].upper() + s[1:] if s else ""

    # Calculate the total energy quantites
    # NOTE: We must generate a new dictionary in order not to get
    #       "dictionary changed size during iteration"
    totEnergies = {}
    for key in energies.keys():
        if "polAvg" in key:
            # Strip polAvg, and cast first letter to lowercase
            quantity = key.replace("polAvg", "")
            totEnergies["tot{}".format(quantity)] =\
                          energies[{"polAvg{}"}.format(quantity)] +\
                          energies[{"fluct{}"} .format(quantity)]
            # Protect data
            totEnergies["tot{}".format(firstCharToUpper(quantity))].\
                          setflags(write=False)

    # Merge dicts
    energies.update(totEnergies)

    # Convert to physical
    if convertToPhysical:
        for key in energies.keys():
            if key[-2:] == "EE":
                varName = "eEnergy"
            elif key[-2:] == "EI":
                varName = "iEnergy"
            elif key == "t":
                varName = "t"

            energies[key] = uc.physicalConversion(energies[key], varName)

    return energies
#}}}

def calcEnergiesWPy():
    # Collect
    varStrings = ("lnN",)

    variables  = collectiveCollect(paths, varStrings)
    lnN        = variables["lnN"]
    n          = calcN(lnN)
    # Delete objects not needed to avoid unnecessary memory consumption
    del variables, lnN

    nPolAvg = polAvg(n)
    nFluct  = n - nPolAvg

    # Delete objects not needed to avoid unnecessary memory consumption
    del n

    polAvgPerpKinEE, polAvgPerpKinEI, fluctPerpKinEE, fluctPerpKinEI =\
            calcPerpEnergiesWPy(n)

    # FIXME: Remember to collect time
    YOU ARE HERE:  Parallel energies

    momDensPar = variables["momDensPar"]
    jPar       = variables["jPar"]

                  ,\
                  "jPar",\

    varStrings = (\
                  "lnN",\
                  "phi",\
                  "momDensPar",\
                  "jPar",\
                 )
    variables  = collectiveCollect(paths, varStrings)
    phi        = variables["phi"]
    # Fluct

                 "t_array"        ,\
    clean with del big arrays





# FIXME: Add docstring (maybe just write it directly into the thesis)
# FIXME: Input nPolAvg and nFluct
#{{{calcPerpEnergiesWPy
def calcPerpEnergiesWPy(n):
    # Collect
    varStrings = ("phi",)
    variables  = collectiveCollect(paths, varStrings)
    phi        = variables["phi"]
    phiPolAvg  = polAvg(phi)
    phiFluct   = phi - phiPolAvg

    # Delete the dictionary and phi to avoid unnecessary memory consumption
    del variables, phi

    uPerpPolAvgSquared =\
            calcGradPerpADotGradPerpB(phiPolAvg, phiPolAvg)

    uPerpFluctUPerpPolAvg =\
            calcGradPerpADotGradPerpB(phiFluct, phiFluct)

    uPerpFluctSquared =\
            calcGradPerpADotGradPerpB(phiFluct, phiFluct)

    nPolAvguPerpFluctSquaredPolAvg = polAvg(nPolAvg*uPerpFluctSquared)
    nFluctuPerpFluctUPerpPolAvgTimesTwoPolAvg =\
            polAvg(nFluctuPerpFluctUPerpPolAvgTimesTwoPolAvg)

    EPerpFluctIntegrand =\
            nFluct*uPerpPolAvgSquared +\
            2.0*nPolAvg*uPerpFluctUPerpPolAvg +\
            (nPolAvg*uPerpFluctSquared - nPolAvguPerpFluctSquaredPolAvg) +\
            (2.0*nFluct*uPerpFluctUPerpPolAvg -\
                    nFluctuPerpFluctUPerpPolAvgTimesTwoPolAvg)

    # Delete variables to avoid unnecessary memory consumption
    del uPerpFluctSquared, uPerpFluctUPerpPolAvg

    # NOTE: This is not the energy per unit volume as we have not taken
    #       the poloidal average of nFluct*uPerpFluctSquared.
    # NOTE: int <a> d theta = int a d theta
    EPerpPolAvgIntegrand =\
            nPolAvg*uPerpPolAvgSquared +\
            nFluct*uPerpFluctSquared +\
            nPolAvguPerpFluctSquaredPolAvg +\
            nFluctuPerpFluctUPerpPolAvgTimesTwoPolAvg

    EPerpFluct  = volInt(EPerpFluctIntegrand)
    EPerpPolAvg = volInt(EPerpPolAvgIntegrand)

    # Delete variables to avoid unnecessary memory consumption
    del EPerpFluctIntegrand, EPerpPolAvgIntegrand

    # NOTE: All collected variables are normalized, so the normalized
    #       energy is the output
# FIXME: Collect mu
    polAvgPerpKinEE = (1.0/mu)*EPerpPolAvg
    polAvgPerpKinEI = EPerpPolAvg
    fluctPerpKinEE  = (1.0/mu)*EPerpFluct
    fluctPerpKinEI  = EPerpFluct

    return polAvgPerpKinEE, polAvgPerpKinEI, fluctPerpKinEE, fluctPerpKinEI
#}}}

# FIXME
invRhoSquared = (1/dh.rho)**(2.0)
#{{{calcGradPerpADotGradPerpB
def calcGradPerpADotGradPerpB(a,b):
    #{{{docstring
    r"""
    Returns \grad_\perp a \cdot \grad_\perp b

    Parameters
    ----------
    a : array-like
        The a in the expression above.
    b : array-like
        The b in the expression above.

    Returns
    -------
    gradPerpADotGradPerpB
    """
    #}}}
    return DDX(a)*DDX(b) + invRhoSquared*DDZ(a)*DDZ(b)
#}}}

# FIXME: Move to collectAndCalcHelpers
#{{{volInt
def volInt(f, dh):
    #{{{docstring
    """
    Returns the volume integrated of a field.

    The integration uses the rectangle rule in 3D, and assumes that the
    point is in the center of a volume cell.
    A simple, stupid, safe approach has been used

    Parameters
    ----------
    f : array 4d
       The variable to volume integrate

    Returns
    -------
    result : array 1d
       The time trace of the volume integrated variable.
    """
    #}}}
    # Initialize the result
    result = 0.0

    _, nx, ny, nz = f.shape
    for xInd in range(nx):
        for yInd in range(ny):
            for zInd in range(nz):
                # NOTE: J = dh.rho
                result +=
                    f       [:, xInd, yInd, zInd] *\
                    dh.rho  [   xInd, yInd]       *\
                    dh.dx   [   xInd, yInd]       *\
                    dh.dy   [   xInd, yInd]       *\
                    dh.dz

    return result
#}}}
