#!/usr/bin/env python

"""
Contains function which calculates the energies
"""

from ..plotHelpers import PlotHelper, collectiveCollect, seqCMap, seqCMap3
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os

#{{{calcEnergy
def calcEnergy(paths                   ,\
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
            tEnd   = tSlice[tCounter].end
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
