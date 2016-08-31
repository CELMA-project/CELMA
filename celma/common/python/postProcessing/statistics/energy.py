#!/usr/bin/env python

"""
Contains function which deals with the post-processing of energies
"""

import matplotlib.pyplot as plt
from .collectiveCollect import collectiveCollect

#{{{collectEnergy
def collectEnergy(paths):
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
        A dictionary of the energies
    """


    varStrings= [\
                 "perpKinEE",\
                 "parKinEE",\
                 "totKinEE",\
                 "perpKinEI",\
                 "parKinEI",\
                 "totKinEI",\
                ]

    energies = collectFromSeveral(paths, varStrings)

    return energies
#}}}

#{{{plotEnergies
def plotEnergies(paths, speciesType, showPlot = False):
    """
    Collects and plots the energies

    Parameters
    ----------
    paths : iterable of strings
        The paths to collect from. Must be in ascending order of the
        simulation time, as the energies are being concatenated
    speciesType : ["electron"|"ion"|"both"]
        What species on should plot for
    showPlot : bool
        Whether or not the plots should be displayed. Default is false

    Return
    ------
    energies : dict
        A dictionary of the energies
    """

    # TODO: Make plots nice

    # Collect the energies
    energies = collectEnergy(path)

    if speciesType == "electron":
        searchStrings = ["EE"]
    elif speciesType == "ion":
        searchStrings = ["EI"]
    elif speciesType == "both":
        searchStrings = ["EE", "EI"]

    for var in energies.keys():
        for searchString in searchStrings:
            if searchString in var:
                plt.plot(energies[var], label=var)
    plt.legend()

    if showPlot:
        plt.show()

    # TODO: Save plot
#}}}
