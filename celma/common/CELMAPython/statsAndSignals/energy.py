#!/usr/bin/env python

"""
Contains function which deals with the post-processing of energies
"""

from ..plotHelpers import PlotHelper, collectiveCollect, seqCMap, seqCMap3
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
                 "perpKinEE"      ,\
                 "parKinEE"       ,\
                 "totKinEE"       ,\
                 "perpKinEI"      ,\
                 "parKinEI"       ,\
                 "totKinEI"       ,\
                 "polAvgPerpKinEE",\
                 "polAvgParKinEE" ,\
                 "polAvgTotKinEE" ,\
                 "polAvgPerpKinEI",\
                 "polAvgParKinEI" ,\
                 "polAvgTotKinEI" ,\
                 "potEE"          ,\
                 "polAvgPotEE"    ,\
                 "totN"           ,\
                 "polAvgTotN"     ,\
                 "t_array"        ,\
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

    #{{{plotKinEnergy
    def plotKinEnergy(self, speciesType):
        #{{{docstring
        """
        Plots the kinetic energies

        For the kinetic energy:
        One figure for each species will be made with the subplots:
            * Perp + par
            * Par
            * Perp
        On each subplot there will be three lines:
            * The full quantity
            * The poloidally averaged quantity
            * The fluctuation quantity.

        Parameters
        ----------
        speciesType : ["electron"|"ion"]
            What species one should plot for
        """
        #}}}

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        gs = GridSpec(nrows=3, ncols=1)
        axes = {"totAx"  : fig.add_subplot(gs[0])}
        axes["parAx" ] = fig.add_subplot(gs[1], sharex=axes["totAx"])
        axes["perpAx"] = fig.add_subplot(gs[2], sharex=axes["totAx"])

        # Get the colors
        colors = seqCMap3(np.linspace(0, 1, 3))

        # Find the keys to plot
        if speciesType == "electrons":
            searchString = "KinEE"
            species = "e"
        elif speciesType == "ions":
            searchString = "KinEI"
            species = "i"
        else:
            message = "speciesType {} not implemented.".format(speciesType)
            raise NotImplementedError(message)

        keys = tuple(key for key in self._energy.keys()\
                if (searchString in key)\
                and ("Units" not in key)\
                and ("Norm" not in key))

        for nr, key in enumerate(keys):
            if "tot" in key.lower():
                ax    = axes["totAx"]
                label = "{} $+$ {}".format(
                    self._genLeg.format(species, r"\mathrm{kin},\perp"),\
                    self._genLeg.format(species, r"\mathrm{kin},\parallel"))
            elif "perp" in key.lower():
                ax    = axes["perpAx"]
                label = self._genLeg.format(species, r"\mathrm{kin},\perp")
            elif "par" in key.lower():
                ax    = axes["parAx"]
                label = self._genLeg.format(species, r"\mathrm{kin},\parallel")

            if "polAvg" in key:
                color = colors[1]
                # Reset label
                labels = label.split("_")
                # Add <>
                # First character is $
                labels[0] = r"{}\langle {} \rangle".\
                        format(labels[0][0], labels[0][1] )
                # Join the labels
                label  = "_".join(labels)
            elif "fluct" in key:
                color = colors[2]
                # Reset label
                labels = label.split("_")
                # Add tilde
                # First character is $
                labels[0] = r"{}\tilde{{{}}}".\
                        format(labels[0][0], labels[0][1] )
                # Join the labels
                label  = "_".join(labels)
            else:
                color = colors[0]

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
                       "{}KinEnergy".format(species)),\
                       self._extension)
            self._helper.savePlot(fig, fileName)

        plt.close(fig)
    #}}}

    #{{{plotPotEnergy
    def plotPotEnergy(self):
        #{{{docstring
        """
        Plots the potential energies

        For the potential energy:
        There is no split in perpendicular and parallel.
        As the ion temperature is zero, all the potential energy will
        come from the electrons:
            * Total ion and electron kinetic vs potential full quantity
            * Total ion and electron kinetic vs potential poloidal
              averaged quantity
            * Total ion and electron kinetic vs potential fluctuation quantity
        """
        #}}}

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        gs = GridSpec(nrows=3, ncols=1)
        axes = {"fullAx"  : fig.add_subplot(gs[0])}
        axes["avgAx"  ] = fig.add_subplot(gs[1], sharex=axes["fullAx"])
        axes["fluctAx"] = fig.add_subplot(gs[2], sharex=axes["fullAx"])

        # Get the colors
        colors = seqCMap(np.linspace(0, 1, 3))

        lines = {}
        for key, ax in axes.items():
            if "full" in key:
                #{{{Full variable
                lines["el" ] = {"line" : self._energy["totKinEE"],\
                                "color": colors[0],\
                                "label":\
                                    (r"$E_{e, \mathrm{kin}, \perp} + "
                                     r" E_{e, \mathrm{kin}, \parallel}$"),\
                                }

                lines["ion"] = {"line" : self._energy["totKinEI"],\
                                "color": colors[1],\
                                "label":\
                                    (r"$E_{i, \mathrm{kin}, \perp} + "
                                     r" E_{i, \mathrm{kin}, \parallel}$"),\
                                }

                lines["pot"] = {"line" : self._energy["potEE"],\
                                "color": colors[2],\
                                "label": r"$E_{e, \mathrm{pot}, \perp}$",\
                                }
                #}}}
            elif "avg" in key:
                #{{{Poloidally averaged variable
                lines["el" ] = {"line" : self._energy["polAvgTotKinEE"],\
                                "color": colors[0],\
                                "label":\
                      (r"$\langle E \rangle_{e, \mathrm{kin}, \perp} + "
                       r" \langle E \rangle_{e, \mathrm{kin}, \parallel}$"),\
                                }

                lines["ion"] = {"line" : self._energy["polAvgTotKinEI"],\
                                "color": colors[1],\
                                "label":\
                      (r"$\langle E \rangle_{i, \mathrm{kin}, \perp} + "
                       r" \langle E \rangle_{i, \mathrm{kin}, \parallel}$"),\
                                }

                lines["pot"] = {"line" : self._energy["polAvgPotEE"],\
                                "color": colors[2],\
                                "label":\
                      r"$\langle E \rangle_{e, \mathrm{pot}, \perp}$",\
                                }
                #}}}
            elif "fluct" in key:
                #{{{Fluctuating variable
                lines["el" ] = {"line" : self._energy["fluctTotKinEE"],\
                                "color": colors[0],\
                                "label":\
                      (r"$\widetilde{E}_{e, \mathrm{kin}, \perp} + "
                       r" \widetilde{E}_{e, \mathrm{kin}, \parallel}$"),\
                                }

                lines["ion"] = {"line" : self._energy["fluctTotKinEI"],\
                                "color": colors[2],\
                                "label":\
                      (r"$\widetilde{E}_{i, \mathrm{kin}, \perp} + "
                       r" \widetilde{E}_{i, \mathrm{kin}, \parallel}$"),\
                                }

                lines["pot"] = {"line" : self._energy["fluctPotEE"],\
                                "color": colors[2],\
                                "label":\
                      r"$\widetilde{E}_{e, \mathrm{pot}, \perp}$",\
                                }
                #}}}

            for _, line in lines.items():
                ax.plot(self._energy["t"]    ,\
                        line["line"]         ,\
                        color = line["color"],\
                        label = line["label"])

        # Turn of x-labels
        axes["fullAx"].tick_params(labelbottom="off")
        axes["avgAx"].tick_params(labelbottom="off")
        # Set axis labels
        axes["avgAx"] .set_ylabel(self._varLabel, labelpad=50)
        axes["fluctAx"].set_xlabel(self._timeLabel)

        # Make the plot look nice
        for key in ["fullAx", "avgAx", "fluctAx"]:
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
                format(os.path.join(self._savePath, "potEnergy"),\
                       self._extension)
            self._helper.savePlot(fig, fileName)

        plt.close(fig)
    #}}}
#}}}
