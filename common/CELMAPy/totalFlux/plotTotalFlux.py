#!/usr/bin/env python

"""Class for totalFlux plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import SizeMaker, plotNumberFormatter, seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotTotalFlux
class PlotTotalFlux(PlotSuperClass):
    """
    Class which contains the total flux data and the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)
    #}}}

    #{{{setData
    def setData(self, totalFluxes, mode, timeAx=True):
        #{{{docstring
        """
        Sets the total fluxes to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
            Dictionary with the keys:
                * "parElIntFlux"  - Time trace of the parallel electron flux
                                    crossing the end plate.
                * "parIonIntFlux" -  Time trace of the parallel ion flux
                                    crossing the end plate.
                * "perpIntFlux"   - Time trace of the perpendicular * flux.
                * "timeIntEl"     - Total number of electrons lost in
                                    the parallel direction.
                * "timeIntIon"    - Total number of ion lost in the parallel
                                    direction.
                * "timeIntPerp"   - Total number of particles lost in the
                                    perpendicular direction.
                * "time"          - Array of the time.
                * "rho"           - The fixed rho value.
                * "z"             - The fixes z value.
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        timeAx : bool
            Whether or not the time should be on the x axis
        """
        #}}}

        self._timeAx = timeAx

        # Set the member data
        self._parElFlux  = totalFluxes.pop("parElIntFlux")
        self._parIonFlux = totalFluxes.pop("parIonIntFlux")
        self._perpFlux   = totalFluxes.pop("perpIntFlux")
        self._parElLoss  = totalFluxes.pop("timeIntEl")
        self._parIonLoss = totalFluxes.pop("timeIntIon")
        self._perpLoss   = totalFluxes.pop("timeIntPerp")
        self._t          = totalFluxes.pop("time")
        self._rho        = totalFluxes.pop("rho")
        self._z          = totalFluxes.pop("z")
        self._mode       = mode

        # Obtain the color (pad away brigthest colors)
        pad = 1
        self._colors = seqCMap3(np.linspace(0, 1, 3+pad))

        self._prepareLabels()

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
            "{}{}".format("totalFluxes", self._fluctName))

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

        # Set label templates
        if self.uc.convertToPhysical:
            self._units = "$[s^{-1}]$"
            normalization = ""
        else:
            self._units = "$[]$"
            normalization = r"/n_0c_s\rho_s^{2}"

        templatePlate  = r"$\iint {}{} \rho\mathrm{{d}}\theta\mathrm{{d}}\rho$"
        templateOuterR = r"$\iint {}{} \rho\mathrm{{d}}\theta\mathrm{{d}}z$"

        if self._mode == "normal":
            var = "nu"
            self._fluctName = ""
        elif self._mode == "fluct":
            var = "\widetilde{n}\widetilde{u}"
            self._fluctName = "-fluct"
        else:
            message = "'{}'-mode not implemented.".format(self._mode)
            raise NotImplementedError(message)

        elTxt   = "{}_{{e,\parallel}}".format(var)
        ionTxt  = "{}_{{i,\parallel}}".format(var)
        perpTxt = "{}_{{E,\perp}}"    .format(var)

        self._elLegend   = templatePlate .format(elTxt  , normalization)
        self._ionLegend  = templatePlate .format(ionTxt , normalization)
        self._perpLegend = templateOuterR.format(perpTxt, normalization)

        # Set total number text
        parElLossTxt = \
            ("$\mathrm{{Total\; \parallel\; electron\;flux}}=${}$\;"
             " \mathrm{{particles}}$").\
                format(plotNumberFormatter(self._parElLoss, None))
        parIonLossTxt = \
            ("$\mathrm{{Total\; \parallel\; ion\; \;flux}}=${}$\;"
             " \mathrm{{particles}}$").\
                format(plotNumberFormatter(self._parIonLoss, None))
        self._parElIonLossTxt = "{}\n{}".format(parElLossTxt, parIonLossTxt)
        self._perpLossTxt = \
            ("$\mathrm{{Total\; \perp\; \;flux}}=2\cdot${}$\;"
             " \mathrm{{particles}}$").\
                format(plotNumberFormatter(self._perpLoss, None))

        # Set values
        self._ph.rhoTxtDict["value"] = plotNumberFormatter(self._rho, None)
        self._ph.zTxtDict  ["value"] = plotNumberFormatter(self._z  , None)

        self._parTitle = r"{}".\
            format(self._ph.zTxtDict["constZTxt"].format(self._ph.zTxtDict))

        self._perpTitle = r"{}".\
            format(self._ph.rhoTxtDict["constRhoTxt"].\
                   format(self._ph.rhoTxtDict))

        # Set the time label
        if self._timeAx:
            self._timeLabel = self._ph.tTxtDict["tTxtLabel"]
        else:
            self._timeLabel = "$\mathrm{Time}$ $\mathrm{index}$"
    #}}}

    #{{{plotSaveShowTotalFlux
    def plotSaveShowTotalFlux(self):
        """
        Performs the actual plotting.

        setData and setDataForSecondVar needs to be called before
        calling this function.
        """

        textPos = (0.96, 0.91)

        # Create the plot
        figSize = SizeMaker.standard(w=5, a=1)
        fig, (elIonAx, perpAx) =\
                plt.subplots(nrows=2, figsize=figSize, sharex=True)

        # Plot the parallel integrated flux
        elIonAx.plot(self._t                ,\
                     self._parElFlux        ,\
                     color = self._colors[0],\
                     label = self._elLegend ,\
                     alpha = 0.7            ,\
                    )
        elIonAx.plot(self._t                ,\
                     self._parIonFlux       ,\
                     color = self._colors[2],\
                     label = self._ionLegend,\
                     alpha = 0.7            ,\
                    )
        # Get the textSize
        txtSize = elIonAx.legend().get_texts()[0].get_fontsize()
        # Add text
        lossTxtElIon = elIonAx.\
                        text(*textPos,\
                            self._parElIonLossTxt,\
                            transform = elIonAx.transAxes,\
                            ha="right", va="top",\
                            bbox={"facecolor":"white", "alpha":0.5, "pad":10},\
                            )
        lossTxtElIon.set_fontsize(txtSize)

        # Plot the perpendicular integrated flux
        perpAx.plot(self._t                 ,\
                    self._perpFlux          ,\
                    color = self._colors[1] ,\
                    label = self._perpLegend,\
                   )

        # Add text
        lossTxtElIon = perpAx.\
                        text(*textPos,\
                            self._perpLossTxt,\
                            transform = perpAx.transAxes,\
                            ha="right", va="top",\
                            bbox={"facecolor":"white", "alpha":0.5, "pad":10},\
                            )
        lossTxtElIon.set_fontsize(txtSize)

        # Set decorations
        elIonAx.set_title(self._parTitle)
        perpAx .set_title(self._perpTitle)

        elIonAx.set_ylabel(self._units)
        perpAx .set_ylabel(self._units)

        perpAx.set_xlabel(self._timeLabel)

        # Make the plot look nice
        self._ph.makePlotPretty(elIonAx, loc="lower right", ybins = 6)
        self._ph.makePlotPretty(perpAx, loc="lower right",\
                                rotation = 45, ybins = 6)

        # Adjust the subplots
        fig.subplots_adjust(hspace=0.2, wspace=0.35)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
