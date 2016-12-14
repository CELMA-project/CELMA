#!/usr/bin/env python

"""
Contains a class used to plot the energy
"""

from ..plotHelpers import PlotHelper, collectiveCollect, seqCMap, seqCMap3
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os

#{{{PlotEnergy
class PlotEnergy(PlotsSuperClass):
    """
    Class which contains the energy data and the plotting configuration.
    """

    #{{{__init___
    def __init__(self    ,\
                 *args   ,\
                 energies,\
                 **kwargs,\
                 ):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets the member data
        3. Prepares the labels

        Parameters
        ----------
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

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the member data
        self._energies = energies

        # Get the units (eventually convert to physical units)
        # NOTE: Need to cast to a tuple to avoid
        #       "RuntimeError: dictionary changed size during iteration"
        for key in tuple(self._energies.keys()):
            # Ions and electrons are normalized in the same manner and
            # have the same units
            norm  = self.uc.conversionDict["eEnergy"]["normalization"]
            units = self.uc.conversionDict["eEnergy"]["units"]

        # Set the variable label
        if self.convertToPhysical:
            self._varLabel = r"$\mathrm{{Energy}}$ $[{}]$".format(units)
        else:
            self._varLabel = r"$\mathrm{{Energy}}{}$".format(norm)

        # Set the time label
        self._timeLabel = self.ph.tTxtDict["tTxtLabel"].\
                          format(self.uc.conversionDict["t"])

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
        axes = {"sumAx"  : fig.add_subplot(gs[0])}
        axes["parAx" ] = fig.add_subplot(gs[1], sharex=axes["sumAx"])
        axes["perpAx"] = fig.add_subplot(gs[2], sharex=axes["sumAx"])

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

        keys = tuple(key for key in self._energies.keys()\
                     if (searchString in key))

        for nr, key in enumerate(keys):
            if "sum" in key.lower():
                ax    = axes["sumAx"]
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

            ax.plot(self._energies["t"],\
                    self._energies[key],\
                    color = color    ,\
                    label = label)

        # Turn of x-labels
        axes["sumAx"].tick_params(labelbottom="off")
        axes["parAx"].tick_params(labelbottom="off")
        # Set axis labels
        axes["parAx"] .set_ylabel(self._varLabel, labelpad=50)
        axes["perpAx"].set_xlabel(self._timeLabel)

        # Make the plot look nice
        for key in ["sumAx", "parAx", "perpAx"]:
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
                lines["el" ] = {"line" : self._energies["sumKinEE"],\
                                "color": colors[0],\
                                "label":\
                                    (r"$E_{e, \mathrm{kin}, \perp} + "
                                     r" E_{e, \mathrm{kin}, \parallel}$"),\
                                }

                lines["ion"] = {"line" : self._energies["sumKinEI"],\
                                "color": colors[1],\
                                "label":\
                                    (r"$E_{i, \mathrm{kin}, \perp} + "
                                     r" E_{i, \mathrm{kin}, \parallel}$"),\
                                }

                lines["pot"] = {"line" : self._energies["potEE"],\
                                "color": colors[2],\
                                "label": r"$E_{e, \mathrm{pot}, \perp}$",\
                                }
                #}}}
            elif "avg" in key:
                #{{{Poloidally averaged variable
                lines["el" ] = {"line" : self._energies["polAvgSumKinEE"],\
                                "color": colors[0],\
                                "label":\
                      (r"$\langle E \rangle_{e, \mathrm{kin}, \perp} + "
                       r" \langle E \rangle_{e, \mathrm{kin}, \parallel}$"),\
                                }

                lines["ion"] = {"line" : self._energies["polAvgSumKinEI"],\
                                "color": colors[1],\
                                "label":\
                      (r"$\langle E \rangle_{i, \mathrm{kin}, \perp} + "
                       r" \langle E \rangle_{i, \mathrm{kin}, \parallel}$"),\
                                }

                lines["pot"] = {"line" : self._energies["polAvgPotEE"],\
                                "color": colors[2],\
                                "label":\
                      r"$\langle E \rangle_{e, \mathrm{pot}, \perp}$",\
                                }
                #}}}
            elif "fluct" in key:
                #{{{Fluctuating variable
                lines["el" ] = {"line" : self._energies["fluctSumKinEE"],\
                                "color": colors[0],\
                                "label":\
                      (r"$\widetilde{E}_{e, \mathrm{kin}, \perp} + "
                       r" \widetilde{E}_{e, \mathrm{kin}, \parallel}$"),\
                                }

                lines["ion"] = {"line" : self._energies["fluctSumKinEI"],\
                                "color": colors[2],\
                                "label":\
                      (r"$\widetilde{E}_{i, \mathrm{kin}, \perp} + "
                       r" \widetilde{E}_{i, \mathrm{kin}, \parallel}$"),\
                                }

                lines["pot"] = {"line" : self._energies["fluctPotEE"],\
                                "color": colors[2],\
                                "label":\
                      r"$\widetilde{E}_{e, \mathrm{pot}, \perp}$",\
                                }
                #}}}

            for _, line in lines.items():
                ax.plot(self._energies["t"]    ,\
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
