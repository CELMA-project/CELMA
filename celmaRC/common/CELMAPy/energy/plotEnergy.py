#!/usr/bin/env python

"""Class for energy plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotEnergy
class PlotEnergy(PlotSuperClass):
    """
    Class which contains the energy data and the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (27,12), **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor

        Parameters
        ----------
        pltSize : tuple
            The size of the plot
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set the plot size
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, energies, timeAx=True):
        #{{{docstring
        """
        Sets the energies to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        energies : dict
            Dictionary with the keys:
            * "perpKinEE" - Time trace of perpendicular kinetic electron energy
            * "parKinEE"  - Time trace of parallel kinetic electron energy
            * "perpKinEI" - Time trace of perpendicular kinetic ion energy
            * "parKinEI"  - Time trace of parallel kinetic ion energy
            * "polEE"     - Time trace of potential electron energy
            * "time"      - Time corresponding time
        timeAx : bool
            Whether or not the time should be on the x axis
        """
        #}}}

        self._timeAx = timeAx

        # Set the member data
        self._energies = energies

        # Obtain the colors. One for the electrons, one for the ions
        # NOTE: Pad away brigthest colors
        pad = 1
        self._colors = seqCMap3(np.linspace(0, 1, 2+pad))

        self._prepareLabels()

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath, "energies")

        if not(timeAx):
            self._fileName += "Indices"
        if (self._sliced):
            self._fileName += "Sliced"

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Set var label template
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
        else:
            unitsOrNormalization = "${normalization}$"

        self._varLabelTemplate = r"$\mathrm{{E}}_{{{},{}{}}}$ " + unitsOrNormalization

        # Set the time label
        if self._timeAx:
            self._timeLabel = self._ph.tTxtDict["tTxtLabel"]
        else:
            self._timeLabel = "$\mathrm{Time}$ $\mathrm{index}$"

        # Make a dict for easy plotting of perp and par
        self._perpParLatex = {"perp":r",\perp", "par":r",\parallel"}
    #}}}

    #{{{plotSaveShowEnergy
    def plotSaveShowEnergy(self):
        """
        Performs the actual plotting.
        """

        # Create the plot
        # Not sharex="col" as it seems like the ticks are lost when setting
        # noneAx axis off
        fig,\
        ((kinEParAx, kinIParAx),\
         (kinEPerpAx, kinIPerpAx),\
         (potEAx, noneAx)) =\
                plt.subplots(ncols=2, nrows=3, figsize=self._pltSize)
        noneAx.set_axis_off()

        axesToTurnOff = (kinEParAx, kinIParAx, kinEPerpAx)
        for ax in axesToTurnOff:
            ax.tick_params(labelbottom="off")

        axDict = {\
                  "kinEParAx" : kinEParAx ,\
                  "kinIParAx" : kinIParAx ,\
                  "kinEPerpAx": kinEPerpAx,\
                  "kinIPerpAx": kinIPerpAx,\
                  "potEAx"    : potEAx    ,\
                  }

        # Plot the kinetic energies
        for eOrI, color in zip(("e","i"), self._colors):
            for perpOrPar in ("perp", "par"):
                axKey  = "kin{}{}Ax".format(eOrI.capitalize(),\
                                        perpOrPar[0].capitalize()+perpOrPar[1:])
                varKey = "{}KinE{}".format(perpOrPar, eOrI.capitalize())

                if self._timeAx:
# FIXME:
                    print("var={}".format(len(self._energies["time"])))
                    print("time={}".format(len(self._energies[varKey])))
# FIXME:
                    axDict[axKey].plot(self._energies["time"],\
                                       self._energies[varKey],\
                                       color=color)
                else:
                    axDict[axKey].plot(self._energies[varKey],\
                                       color=color)

                # Get the y-label
                varLabel = self._varLabelTemplate.format(\
                                        r"\mathrm{kin}",\
                                        eOrI,\
                                        self._perpParLatex[perpOrPar],\
                                        **self.uc.conversionDict[eOrI+"Energy"])
                axDict[axKey].set_ylabel(varLabel)

        # Plot the potential energies
        if self._timeAx:
            potEAx.plot(self._energies["time"],\
                        self._energies["potEE"],\
                        color=self._colors[0])
        else:
            potEAx.plot(self._energies["potEE"],\
                        color=self._colors[0])
        varLabel = self._varLabelTemplate.format(\
                                r"\mathrm{pot}", "e", "",\
                                **self.uc.conversionDict["eEnergy"])
        potEAx.set_ylabel(varLabel)

        # Set x labels
        potEAx    .set_xlabel(self._timeLabel)
        kinIPerpAx.set_xlabel(self._timeLabel)

        # Set titles
        kinEParAx.set_title("$\mathrm{Electron \quad energy}$")
        kinIParAx.set_title("$\mathrm{Ion \quad energy}$")

        # Make the plot look nice
        for key in axDict.keys():
            # Make the ax pretty
            self._ph.makePlotPretty(axDict[key],\
                                    yprune   ="both",\
                                    rotation = 45   ,\
                                    ybins    = 6    ,\
                                    legend   = False,\
                                    )

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.1, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
