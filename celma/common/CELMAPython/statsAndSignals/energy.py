#!/usr/bin/env python

"""
Contains function which deals with the post-processing of energies
"""

from ..plotHelpers import PlotHelper, collectiveCollect, seqCMap3
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os

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
        A dictionary of the energies (including the time)
    """


    varStrings= (\
                 "perpKinEE",\
                 "parKinEE",\
                 "totKinEE",\
                 "perpKinEI",\
                 "parKinEI",\
                 "totKinEI",\
                 "t_array",\
                )

    energies = collectiveCollect(paths, varStrings)

    # Change keyname of t_array to t
    energies["t"] = energies.pop("t_array")

    return energies
#}}}

#{{{PlotEnergy
class PlotEnergy(object):
    """
    Class which contains the energy data and the plotting configuration.
    """

    #{{{__init___
    def __init__(self                     ,\
                 paths                    ,\
                 energy                   ,\
                 convertToPhysical = False,\
                 showPlot          = False,\
                 savePlot          = False,\
                 extension         = "png",\
                 savePath          = "."  ,\
                 pltSize           = None ,\
                 ):
        #{{{docstring
        """
        The constructor for the PlotEnergy object.

        Sets the member data.

        Parameters
        ----------
        paths : str
            Paths to collect from (used to make the PlotHelper object)
        energy : dictionary
            Contains the energy in the following keys:
                * perpKinEE - The perpendicular kinetic electron energy
                * parKinEE  - The paralell kinetic electron energy
                * totKinEE  - The total kinetic electron energy
                * perpKinEI - The perpendicular kinetic ion energy
                * parKinEI  - The parallel kinetic ion energy
                * totKinEI  - The total kinetic ion energy
                * t_array   - The time
        showPlot : bool
            If the plots should be displayed.
        savePlot : bool
            If the plots should be saved.
        extension : str
            Extension to use on the plots
        savePath : str
            Path to save destination. Must exist.
        pltSize : tuple
            Size of the plots given as (x, y)
        """
        #}}}

        # Set the member data
        self._energy    = energy
        self._showPlot  = showPlot
        self._savePlot  = savePlot
        self._extension = extension
        self._savePath  = savePath

        # Get the colors
        self._colors = seqCMap3(np.linspace(0, 1, 3))

        # Make the PlotHelper object
        self._helper = PlotHelper(paths[0]                              ,\
                                  # Copy the array as we do not want to
                                  # share memory
                                  t                 = energy["t"].copy(),\
                                  xguards           = False             ,\
                                  yguards           = False             ,\
                                  convertToPhysical = convertToPhysical ,\
                                 )

        # Get the units (eventually convert to physical units)
        # NOTE: Need to cast to a tuple to avoid
        #       "RuntimeError: dictionary changed size during iteration"
        for key in tuple(self._energy.keys()):
            self._energy[key],\
            self._energy[key+"Norm"],\
            self._energy[key+"Units"] =\
                self._helper.physicalUnitsConverter(self._energy[key], key)

        # Set the variable label
        if self._helper.convertToPhysical:
            self._varLabel = r"$\mathrm{{Energy}}$ $[{}]$".\
                                  format(self._energy["totKinEIUnits"])
        else:
            self._varLabel = r"$\mathrm{{Energy}}{}$".\
                                  format(self._energy["totKinEINorm"])

        # Set the time label
        self._timeLabel = self._helper.tTxtDict["tTxtLabel"].\
                          format(self._helper.tTxtDict)

        # Set the generic legend
        self._genLeg = r"$E_{{{}, {}}}$"
        # Set the plot size
        self._pltSize = pltSize
    #}}}

    #{{{plotEnergies
    def plotEnergies(self, speciesType):
        """
        Plots the energies

        Parameters
        ----------
        speciesType : ["electron"|"ion"]
            What species one should plot for
        """

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        gs = GridSpec(nrows=3, ncols=1)
        axes = {"totAx"  : fig.add_subplot(gs[0])}
        axes["parAx" ] = fig.add_subplot(gs[1], sharex=axes["totAx"])
        axes["perpAx"] = fig.add_subplot(gs[2], sharex=axes["totAx"])

        if speciesType == "electrons":
            searchString = "EE"
            species = "e"
        elif speciesType == "ions":
            searchString = "EI"
            species = "i"
        else:
            message = "speciesType {} not implemented.".format(speciesType)
            raise NotImplementedError(message)

        # Find the keys to plot
        keys = tuple(key for key in self._energy.keys()\
                if (searchString in key)\
                and ("Units" not in key)\
                and ("Norm" not in key))

        for nr, key in enumerate(keys):
            if "tot" in key:
                ax    = axes["totAx"]
                label = self._genLeg.format(species, r"\mathrm{kin, tot}")
                color = self._colors[0]
            elif "par" in key:
                ax    = axes["parAx"]
                label = self._genLeg.format(species, r"\mathrm{kin},\parallel")
                color = self._colors[1]
            elif "perp" in key:
                ax    = axes["perpAx"]
                label = self._genLeg.format(species, r"\mathrm{kin},\perp")
                color = self._colors[2]
            ax.plot(self._energy["t"],\
                    self._energy[key],\
                    color = color    ,\
                    label = label)

        # Turn of x-labels
        axes["totAx"].tick_params(labelbottom="off")
        axes["parAx"].tick_params(labelbottom="off")
        # Set axis labels
        axes["parAx"] .set_ylabel(self._varLabel, labelpad=50)
        axes["perpAx"].set_xlabel(self._timeLabel)

        # Make the plot look nice
        for key in ["totAx", "parAx", "perpAx"]:
            self._helper.makePlotPretty(axes[key]              ,\
                                        yprune   = "both"      ,\
                                        rotation = 45          ,\
                                        loc      = "lower right")

        # Adjust the subplots
        fig.subplots_adjust(hspace=0)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            fileName = "{}.{}".\
                format(os.path.join(self._savePath,\
                       "{}Energy".format(species)),\
                       self._extension)
            self._helper.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
