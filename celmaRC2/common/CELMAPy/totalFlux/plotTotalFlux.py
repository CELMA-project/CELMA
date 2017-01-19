#!/usr/bin/env python

"""Class for totalFlux plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotTotalFlux
class PlotTotalFlux(PlotSuperClass):
    """
    Class which contains the total flux data and the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (15,10), **kwargs):
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
    def setData(self, totalFluxes, mode, timeAx=True):
        #{{{docstring
        """
        Sets the total fluxes to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        totalFluxes : dict
            Dictionary on the form:
            {"parElIntFlux"  : parallelElectronIntegratedFlux,
             "parIonIntFlux" : parallelIonIntegratedFlux,
             "perpIntFlux"   : perpendicularIntegratedFlux,
             "time"          : time,
             "rho"           : rho,
             "z"             : z,
            }
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

        # Set var label template
        if self.uc.convertToPhysical:
            self._units = "$[s^{-1}]$"
            normalization = ""
        else:
            self._units = "$[]$"
            normalization = r"/n_0c_s\rho_s^{2}"

        template = r"$\iint {}{} \rho\mathrm{{d}}\theta\mathrm{{d}}z$"

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

        self._elLegend   = template.format(elTxt  , normalization)
        self._ionLegend  = template.format(ionTxt , normalization)
        self._perpLegend = template.format(perpTxt, normalization)

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

        # Create the plot
        fig, (elIonAx, perpAx) =\
                plt.subplots(nrows=2, figsize=self._pltSize, sharex=True)

        # Plot the parallel integrated flux
        elIonAx.plot(self._t                ,\
                     self._parElFlux        ,\
                     color = self._colors[0],\
                     label = self._elLegend ,\
                    )
        elIonAx.plot(self._t                ,\
                     self._parIonFlux       ,\
                     color = self._colors[2],\
                     label = self._ionLegend,\
                    )

        # Plot the perpendicular integrated flux
        perpAx.plot(self._t                 ,\
                    self._perpFlux          ,\
                    color = self._colors[1] ,\
                    label = self._perpLegend,\
                   )

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
