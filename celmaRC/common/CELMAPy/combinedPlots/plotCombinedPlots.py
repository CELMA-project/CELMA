#!/usr/bin/env python

"""Class for combinedPlots plot"""

from ..superClasses import PlotSuperClass
from ..plotHelpers import plotNumberFormatter, seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotCombinedPlots
class PlotCombinedPlots(PlotSuperClass):
    """
    Class which contains the combined plots data and the plotting configuration.
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
    def setData(self, tt, rf, PDF, PSD, mode):
        #{{{docstring
        """
        Sets the combined plots to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        tt : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:timeTrace, "time":time}
        radialFluxes : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:radialFlux, "time":time}
        PDF : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"pdfX":pdfX, "pdfY":"pdfY"}
        PSD : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {"pdfX":pdfX, "pdfY":"pdfY"}
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        """
        #}}}

        # Set the member data
        self._tt   = tt
        self._rf   = rf
        self._PDF  = PDF
        self._PSD  = PSD
        self._mode = mode

        # Obtain the varname
        ind  = tuple(tt.keys())[0]
        keys = tt[ind].keys()
        self._varName = tuple(var for var in keys if var != "time")[0]
        self._radialFluxName = tuple(var for var in keys if var != "time")[0]
        self._pltVarName = self._ph.getVarPltName(self._varName)


        # Make the keys from the time trace
        self._keys = sorted(self._tt.keys(), reversed = True)

        # To be used in the plot name
        if self._mode == "normal":
            self._fluctName = ""
        elif self._mode == "fluct":
            self._fluctName = "fluct"
        else:
            message = "'{}'-mode not implemented.".format(self._mode)
            raise NotImplementedError(message)

        # Obtain the color (pad away brigthest colors)
        pad = 3
        self._colors = seqCMap3(np.linspace(0, 1, len(combinedPlots.keys())+pad))

        # Prepare the labels
        self._prepareLabels()

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, "combinedPlots", self._fluctName))

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Make vertical titles
        self._verticalTitles = self._makeVerticalTitles()

        # Set var label template
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
        else:
            unitsOrNormalization = "${normalization}$"

        # Make the horizontal titles
        self._horizontalTitles =\
            self._makeHorizontalTitles(unitsOrNormalization)

        # Make the xLabels
        self._xLabels = self._makeXLabels(unitsOrNormalization)
    #}}}

    #{{{_makeVerticalTitles
    def _makeVerticalTitles(self):
        #{{{docstring
        """
        Makes the vertical titles

        Returns
        -------
        verticalTitles : tuple
            Title with the vertical titles (i.e. the positions) in the
            order in descending order.
        """
        #}}}
        verticalTitles = []
        for key, color in zip(self._keys, self._colors):
            # Make the label
            rho, theta, z = key.split(",")

            # Set values
            self._ph.rhoTxtDict  ["value"] =\
                    plotNumberFormatter(float(rho), None)
            self._ph.zTxtDict    ["value"] =\
                    plotNumberFormatter(float(z), None)
            self._ph.thetaTxtDict["value"] =\
                    plotNumberFormatter(float(theta), None)

            # Make the const values
            verticalTitles.append((r"{}$,$ {}$,$ {}").\
                    format(\
                        self._ph.rhoTxtDict  ["constRhoTxt"].\
                            format(self._ph.rhoTxtDict),\
                        self._ph.thetaTxtDict["constThetaTxt"].\
                            format(self._ph.thetaTxtDict),\
                        self._ph.zTxtDict["constZTxt"].\
                            format(self._ph.zTxtDict),\
                          ))

        return tuple(verticalTitles)
    #}}}

    #{{{_makeHorizontalTitles
    def _makeHorizontalTitles(self, unitsOrNormalization):
        #{{{docstring
        """
        Makes the horizontal titles.

        Parameters
        ----------
        unitsOrNormalization : str
            String containing either the units or the normalization

        Returns
        -------
        verticalTitles : tuple
            Title with the horizontal titles in the order
                0. Time trace
                1. Radial flux
                2. PDF
                3. PSD
        """
        #}}}

        horizontalTitles = []

        ttFirstLine  = r"\mathrm{Time \quad trace}" +"\n"
        rfFirstLine  = r"\mathrm{Raidal \quad flux}"+"\n"
        PDFFirstLine = r"\mathrm{PDF}"+"\n"
        PSDFirstLine = r"\mathrm{PSD}"+"\n"

        # The time trace
        if self._mode == "normal":
            ttSecondLine = r"${{}}${}".format(unitsOrNormalization)
        elif self._mode == "fluct":
            ttSecondLine = r"$\widetilde{{{{{{}}}}}}${}".\
                    format(unitsOrNormalization)

        # The radial flux
        if self._mode == "normal":
            rfSecondLine = r"${{}}u_{{{{E\times B, \rho}}}}${}".\
                                format(unitsOrNormalization)
        elif self._mode == "fluct":
            rfSecondLine = (\
                    r"$\widetilde{{{{{{}}}}}} "
                    r"\widetilde{{{{ u }}}}_{{{{E\times B, \rho}}}}${}").\
                    format(unitsOrNormalization)

        # The probability density function
        PDFSecondLine = self._preparePDFLabel()

        # The power spectral density
        PSDSecondLine = self._preparePSDLabel()

        ttSecondLine  = ttSecondLine .format(self._pltvarName)
        rfSecondLine  = rfSecondLine .format(self._pltvarName)
        PDFSecondLine = PDFSecondLine.format(self._pltvarName)
        PSDSecondLine = PSDSecondLine.format(self._pltvarName)

        horizontalTitles.append(ttFirstLine+ttSecondLine)
        horizontalTitles.append(rfFirstLine+rfSecondLine)
        horizontalTitles.append(PDFFirstLine+PDFSecondLine)
        horizontalTitles.append(PSDFirstLine+PSDSecondLine)

        return tuple(horizontalTitles)
    #}}}

    #{{{_preparePDFLabel
    def _preparePDFLabel(self):
        """
        Returns the prepared PDF label
        """

        pdfStr = r"\mathrm{{PDF}}("
        # Get normalOrFluct
        if self._mode == "normal":
            normalOrFluct = r"{{}}"
        elif self._mode == "fluct":
            normalOrFluct = r"\widetilde{{{}}}"
        #  Normalized and non-normalized labels are treated differently
        #  due to the paranthesis and the normalization
        if self.uc.convertToPhysical:
            unitsOrNormalization       = "[1/{units}]"
            preparedPDFLabel = r"$"+\
                               pdfStr +\
                               normalOrFluct +\
                               r")$" +\
                               r" $" + unitsOrNormalization + r"$"
        else:
            unitsOrNormalization       = "{normalization}"
            preparedPDFLabel = r"$"+\
                               pdfStr +\
                               normalOrFluct +\
                               unitsOrNormalization +\
                               r")$"

        return preparedPDFLabel
    #}}}

    #{{{_preparePSDLabel
    def _preparePSDLabel(self):
        """
        Returns the prepared PSD label
        """
        psdStr = r"\mathrm{{PSD}}("

        # Get normalOrFluct
        if self._mode == "normal":
            normalOrFluct = r"{{}}"
        elif self._mode == "fluct":
            normalOrFluct = r"\widetilde{{{}}}"

        #  Normalized and non-normalized labels are treated differently
        #  due to the paranthesis and the normalization
        units = "[\mathrm{{dB}}]"
        if self.uc.convertToPhysical:
            preparedPSDLabel = r"$"+\
                               psdStr +\
                               normalOrFluct +\
                               r")$" +\
                               r" $" + units + r"$"
        else:
            preparedPSDLabel = r"$"+\
                               psdStr +\
                               normalOrFluct +\
                               "{normalization}" +\
                               r")$" +\
                               r" $" + units + r"$"

        return preparedPSDLabel
    #}}}

    #{{{_makeXLabels
    def _makeXLabels(self, unitsOrNormalization):
        #{{{docstring
        """
        Makes the x labels.

        Parameters
        ----------
        unitsOrNormalization : str
            String containing either the units or the normalization.

        Returns
        -------
        verticalTitles : tuple
            Title with the horizontal titles in the order
                0. Time trace
                1. Radial flux
                2. PDF
                3. PSD
        """
        #}}}

        xLabels = []

        timeLabel = self._ph.tTxtDict["tTxtLabel"]

        ttXLabel  = timeLabel
        rfXLabel  = timeLabel
        PSDXLabel = r"$1/t$ $[Hz]$" if self.uc.convertToPhysical\
                    else r"$1/t \omega_{ci}$"

        # Make the PDFXLabel
        PDFXLabelunitsOrNormalization = "[{units}]"
        PDFXLabelunitsOrNormalization = "{normalization}"
        PDFXLabel = r"$" + normalOrFluct + r"$" +\
                    r" $" + xLabelunitsOrNormalization + r"$"
        PDFXLabel.format(unitsOrNormalization)

        xLabels.append(ttXLabel )
        xLabels.append(rfXLabel )
        xLabels.append(PDFXLabel)
        xLabels.append(PSDXLabel)

        return tuple(xLabels)
    #}}}

    @staticmethod
    #{{{turnOffXTicks
    def turnOffXTicks(axes):
        #{{{docstring
        """
        Turns off x-ticks on the input axes

        Parameters
        ----------
        axes : tuple
            Tuple of the axes to turn off the x-ticks
        """
        #}}}
        for ax in axes:
            plt.setp(ax.get_xticklabels(), visible=False)
    #}}}

    #{{{_prepareAxes
    def _prepareAxes(self):
        # http://stackoverflow.com/questions/21191519/python-matplotlib-moving-graph-title-to-the-y-axis
        # http://matplotlib.org/users/gridspec.html

        self._ttAxes, self._rfAxes, self._PDFAxes, self._PSDAxes =\
                self._makeAxes()

        self.turnOffXTicks((*self._ttAxes [:-1],\
                            *self._rfAxes [:-1],\
                            *self._PDFAxes[:-1],\
                            *self._PSDAxes[:-1],\
                           ))

        self._setTitlesAndLables()
    #}}}

    #{{{_makeAxes
    def _makeAxes(self):
        #{{{docsting
        """
        Makes the axes

        Returns
        -------
        ttAxes : axis
            The time trace axis.
        rfAxes : axis
            The radial flux axis.
        PDFAxes: axis
            The probability density function axis.
        PSDAxes: axis
            The power spectral density axis.
        """
        #}}}
        fig = plt.figure()

        gs = gridspec.GridSpec(3, 4,\
                               width_ratios =(3,3,2,2),\
                               height_ratios=(2,2,2,2),\
                               )

        ax1  = plt.subplot(gs[0])
        ax2  = plt.subplot(gs[1])
        ax3  = plt.subplot(gs[2])
        ax4  = plt.subplot(gs[3])
        ax5  = plt.subplot(gs[4] , sharex=ax1, sharey=ax1)
        ax6  = plt.subplot(gs[5] , sharex=ax2, sharey=ax2)
        ax7  = plt.subplot(gs[6] , sharex=ax3, sharey=ax3)
        ax8  = plt.subplot(gs[7] , sharex=ax4, sharey=ax4)
        ax9  = plt.subplot(gs[8] , sharex=ax1, sharey=ax1)
        ax10 = plt.subplot(gs[9] , sharex=ax2, sharey=ax2)
        ax11 = plt.subplot(gs[10], sharex=ax3, sharey=ax3)
        ax12 = plt.subplot(gs[11], sharex=ax4, sharey=ax4)

        ttAxes  = (ax1,ax5,ax9)
        rfAxes  = (ax2,ax6,ax10)
        PDFAxes = (ax3,ax7,ax11)
        PSDAxes = (ax4,ax8,ax12)

        return ttAxes, rfAxes, PDFAxes, PSDAxes,
    #}}}

    #{{{_setTitlesAndLabels
    def _setTitlesAndLabels(self):
        """
        Set the titles and the labels
        """

        # Set horizontal titles
        self._ttAxes[0] .set_title(self._horizontalTitles[0])
        self._rfAxes[0] .set_title(self._horizontalTitles[1])
        self._PDFAxes[0].set_title(self._horizontalTitles[2])
        self._PSDAxes[0].set_title(self._horizontalTitles[3])

        # Set time trace labels
        labelTop    = self._ttAxes[0].set_ylabel("")
        labelMid    = self._ttAxes[1].set_ylabel("")
        labelBottom = self._ttAxes[2].set_ylabel("")

        # Set vertical titles
        titleMid    = self._ttAxes[1].set_title(self._verticalTitles[1])
        titleBottom = self._ttAxes[2].set_title(self._verticalTitles[2])
        # Special case for first axis as ax title is already in use
        titleTop = self._ttAxes[0].\
                    text(0,0  ,\
                        self._verticalTitles[0],\
                        transform=ax1.transAxes, ha="center", va="center",\
                        )
        titleTop.set_fontsize(titleMid.get_fontsize())

        # Magic numbers
        offset = np.array([-0.3, 0.0])

        titleTop    .set_rotation(90)
        titleTop    .set_position(labelTop  .get_position() + offset)
        titleMid    .set_rotation(90)
        titleMid    .set_position(labelMid  .get_position() + offset)
        titleBottom.set_rotation(90)
        titleBottom.set_position(labelBottom.get_position() + offset)
    #}}}

    #{{{plotSaveShowCombinedPlots
    def plotSaveShowCombinedPlots(self):
        """
        Performs the actual plotting.
        """
        for nr, key in enumerate(self._keys):
            self._ttAxes[nr].plot(self._tt[key]["time"]       ,\
                                  self._tt[key][self._varName],\
                                  color=self._colors[nr])

            self._rfAxes[nr].plot(self._rf[key]["time"]              ,\
                                  self._rf[key][self._radialFluxName],\
                                  color=self._colors[nr])

            self._PDFAxes[nr].plot(self._PDF[key]["{}PDFX".format(self._varName)],\
                                   self._PDF[key]["{}PDFY".format(self._varName)],\
                                   color=self._colors[nr])

            # PSD plot
            # Clip the very first point as this is very low
            PSDXVals = self._PSD[key]["{}PSDX".format(self._varName)][1:]
            PSDYVals = np.log10(       self._PSD[key]["{}PSDY".\
                                         format(self._varName)][1:]/\
                                 np.max(self._PSD[key]["{}PSDY".\
                                        format(self._varName)][1:]))
            self._PSDAxes[nr].plot(xVals, yVals, color=color)

        # Use logarithmic scale
        for ax in self._PDFAxes:
            ax.set_yscale("log")
        for ax in self._PSDAxes:
            ax.set_xscale("log")

        # Make the plot look nice
        bottomAxes =\
            (*self._ttAxes [:-1],\
             *self._rfAxes [:-1],\
             *self._PDFAxes[:-1],\
             *self._PSDAxes[:-1])
        topAxes =\
            (*self._ttAxes [0],\
             *self._rfAxes [0],\
             *self._PDFAxes[0],\
             *self._PSDAxes[0])
        for ax in topAxes:
            self._ph.makePlotPretty(ax            ,\
                                    ybins = 3     ,\
                                    legend = False,\
                                    )
        for ax in bottomAxes:
            self._ph.makePlotPretty(ax            ,\
                                    ybins = 3     ,\
                                    rotation = 45 ,\
                                    legend = False,\
                                    )

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName)

        plt.close(fig)
    #}}}
#}}}
