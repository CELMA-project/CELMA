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
def calcEnergy(paths):
    """
    Collects and concatenate the energies

    Parameters
    ----------
    paths : iterable of strings
        The paths to collect from. Must be in ascending order of the
        simulation time, as the variables are being concatenated

    Return
    ------
    energies : dict
        A dictionary of the energies (including the time)
    """


    varStrings= (\
                 "perpKinEE"           ,\
                 "parKinEE"            ,\
                 "sumKinEE"            ,\
                 "perpKinEI"           ,\
                 "parKinEI"            ,\
                 "sumKinEI"            ,\
                 "polAvgPerpKinEE"     ,\
                 "polAvgParKinEE"      ,\
                 "polAvgSumKinEE"      ,\
                 "polAvgPerpKinEI"     ,\
                 "polAvgParKinEI"      ,\
                 "polAvgSumKinEI"      ,\
                 "potEE"               ,\
                 "polAvgPotEE"         ,\
                 "particleNumber"      ,\
                 "t_array"             ,\
                )

    energies = collectiveCollect(paths, varStrings)

    # Change keyname of t_array to t
    energies["t"] = energies.pop("t_array")

    firstCharToLower = lambda s: s[:1].lower() + s[1:] if s else ""
    firstCharToUpper = lambda s: s[:1].upper() + s[1:] if s else ""

    # Calculate fluctuating quantities
    # NOTE: We must generate a new dictionary in order not to get
    #       "dictionary changed size during iteration"
    fluctEnergies = {}
    for key in energies.keys():
        if "polAvg" in key:
            # Strip polAvg, and cast first letter to lowercase
            quantity = firstCharToLower(key.replace("polAvg", ""))
            fluctEnergies["fluct{}".format(firstCharToUpper(quantity))] =\
                          energies[quantity] - energies[key]
            # Protect data
            fluctEnergies["fluct{}".format(firstCharToUpper(quantity))].\
                          setflags(write=False)

    # Merge dicts
    energies.update(fluctEnergies)

    return energies
#}}}
