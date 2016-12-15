#!/usr/bin/env python

""" Contains the PlotHelper class """

from .plotNumberFormatter import plotNumberFormatter
from matplotlib.ticker import MaxNLocator, FuncFormatter
import os

#{{{PlotHelper
class PlotHelper(object):
    """Contains common routines used when making plots"""

    #{{{Static members
    _varPltName = {\
        "lnN"        :  "\ln(n)"          ,\
        "n"          :  "n"               ,\
        "vort"       :  "\Om"             ,\
        "vortD"      :  "\Om^D"           ,\
        "phi"        :  "\phi"            ,\
        "jPar"       :  "j_{\parallel}"   ,\
        "momDensPar" :  "nu_{i,\parallel}",\
        "uIPar"      :  "u_{i,\parallel}" ,\
        "uEPar"      :  "u_{e,\parallel}" ,\
        }
    #}}}

    #{{{__init__
    def  __init__(self                     ,\
                  convertToPhysical = False):
        #{{{docstring
        """
        The constructor for PlotHelper, which:

        * Sets the member data

        Parameters
        ----------
        convertToPhysical : bool
            Whether or not to convert to physical units.
        """
        #}}}

        # Set the member data
        self._convertToPhysical = convertToPhysical
    #}}}

    #{{{makeDimensionStringsDicts
    def makeDimensionStringsDicts(self, unitsConverter):
        #{{{docstring
        """
        Makes the dimension strings dicts.

        The tTxtDict, rhoTxtDict, thetaTxtDict and zTxtDict will be accessable
        from the PlotHelper object, and contain the following keys:
            * units         - String of the units for the variabls
            * normalization - String of the normalization for the variable
            * @Txt          - Variable with possible normalization
            * @TxtLabel     - Label of the variable
            * const@Txt     - Label when the variable is constant in a plot

        Expands the TxTDicts, so that they can easily be used when
        formatting text for labels and titles.
        """
        #}}}

        # String formatting
        self.tTxtDict     =\
            {"normalization":unitsConverter.getNormalization("t"  ),\
             "units"        :unitsConverter.getUnits        ("t"  ) }
        self.rhoTxtDict   =\
            {"normalization":unitsConverter.getNormalization("rho"),\
             "units"        :unitsConverter.getUnits        ("rho") }
        self.zTxtDict     =\
            {"normalization":unitsConverter.getNormalization("z"  ),\
             "units"        :unitsConverter.getUnits        ("z"  ) }
        self.thetaTxtDict = {}

        # Set generic string templates
        self.rhoTxtDict["rhoTxt"] =\
                r"$\rho{0[normalization]}$".format(self.rhoTxtDict)
        self.zTxtDict["zTxt"] =\
                r"$z{0[normalization]}$".format(self.zTxtDict)
        self.thetaTxtDict["constThetaTxt"] =\
                r"$\theta={:d}^{{\circ}}$"
        # Set label and title templates
        if self._convertToPhysical:
            self.rhoTxtDict["rhoTxtLabel"] = "{0[rhoTxt]} $[{0[units]}]$".\
                        format(self.rhoTxtDict)
            self.rhoTxtDict["constRhoTxt"] =\
                        r"{0[rhoTxt]} $=$ {0[value]} ${0[units]}$"
            self.zTxtDict["zTxtLabel"] = "{0[zTxt]} $[{0[units]}]$".\
                        format(self.zTxtDict)
            self.zTxtDict["constZTxt"] =\
                    r"{0[zTxt]} $=$ {0[value]} ${0[units]}$"
            self.tTxtDict["tTxt"] =\
                r"$\mathrm{{t}}{0[normalization]}$ $=$ {0[value]} ${0[units]}$"
            self.tTxtDict["tTxtLabel"] =\
                r"$\mathrm{{t}}{0[normalization]}$ $[{0[units]}]$"
        else:
            self.rhoTxtDict["rhoTxtLabel"] = "{0[rhoTxt]}".\
                                             format(self.rhoTxtDict)
            self.rhoTxtDict["constRhoTxt"] = r"{0[rhoTxt]} $=$ {0[value]}"
            self.zTxtDict  ["zTxtLabel"] = "{0[zTxt]}".\
                                           format(self.zTxtDict)
            self.zTxtDict  ["constZTxt"] = r"{0[zTxt]} $=$ {0[value]}"
            self.tTxtDict["tTxt"] =\
                r"$\mathrm{{t}}{0[normalization]}$ $=$ {0[value]}"
            self.tTxtDict  ["tTxtLabel"] = r"$t{0[normalization]}$"
    #}}}

    @staticmethod
    #{{{getVarPltName
    def getVarPltName(var):
        """
        Routine that returns the variable plot name.
        The return value does not include the $

        Parameters
        ----------
        var : str
            The variable to find the plot name for.

        Retruns
        -------
        varPltName : str
            The plot name
        """

        return PlotHelper._varPltName[var]
    #}}}

    @staticmethod
    #{{{makePlotPretty
    def makePlotPretty(ax                     ,\
                       xprune   = "lower"     ,\
                       yprune   = None        ,\
                       rotation = "horizontal",\
                       loc      = "best"      ,\
                       ):
        """
        Routine that fixes some beauty-mistakes in matplotlib

        Parameters
        ----------
        ax : axis
            The axis to fix.
        xprune : str
            What ticks should be pruned on the x axis.
        yprune : str
            What ticks should be pruned on the y axis.
        rotation : [str | int]
            Rotation of the x axis.
        loc : str
            Location of the legend
        """

        # Avoid silly top value (only for non-log axes)
        try:
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
        except:
            pass
        # Format the tick labels
        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=rotation)
        ax.get_xaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
        ax.get_yaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
        # Plot the legend
        try:
            leg = ax.legend(loc       = loc ,\
                            fancybox  = True,\
                            numpoints = 1   ,\
                            )
            leg.get_frame().set_alpha(0.5)
        except AttributeError as ae:
            if "NoneType" in ae.args[0] and "get_frame" in ae.args[0]:
                print("{0}{1}WARNING: Could not set legend{1}{0}".format("\n","!"*4))
            else:
                raise ae
        # Plot the grid
        ax.grid()
        # Make sure no collision between the ticks
        if ax.get_xscale() != "log":
            ax.xaxis.set_major_locator(MaxNLocator(prune=xprune))

        if ax.get_yscale() != "log":
            # This destroys the ticks on log plots
            ax.yaxis.set_major_locator(MaxNLocator(prune=yprune))
    #}}}

    @staticmethod
    #{{{savePlot
    def savePlot(fig, fileName, extraArtists=None):
        """
        Saves the figure

        Parameters
        ----------
        fig: figure
            The figure.
        fileName : str
            Full path of the plot.
        extraArtist : tuple
            Tuple of bbox_extra_artists to be saved
        """

        # Create path if not exists
        directory = os.path.dirname(fileName)
        if directory != "" and directory != ".":
            if not os.path.exists(directory):
                    os.makedirs(directory)
                    print("{} created".format(directory))


        fig.savefig(fileName,\
                    transparent = True             ,\
                    bbox_inches = "tight"          ,\
                    bbox_extra_artists=extraArtists,\
                    pad_inches  = 0                ,\
                    )

        print("Saved to {}".format(fileName))
    #}}}
#}}}
