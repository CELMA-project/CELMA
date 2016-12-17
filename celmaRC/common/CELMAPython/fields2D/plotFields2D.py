#!/usr/bin/env python

"""
Contains functions for plotting the 2D fields
"""

import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..superClasses import Plot2DSuperClass
from ..plotHelpers import plotNumberFormatter

#{{{Plot2DPerp
class Plot2DPerp(Plot2DSuperClass):
    """
    Class for 2D perp plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (10,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot2DPerp

        * Creates the figure and axes

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        pltSize : tuple
            The size of the plot
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Create figure and axes
        pltSize = (20,15)
        self._fig = plt.figure(figsize = pltSize)
        self._perpAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._perpAx).\
                 append_axes('right', '5%', '5%')
        self._perpAx.grid(True)
    #}}}

    #{{{setData
    def setData(self, X_RT, Y_RT, Z_RT, time, constZ):
        #{{{docstring
        """
        Sets the data to  be plotted

        Parameters
        ----------
        X_RT : array
            A 2d mesh of the Cartesian x coordinates.
        Y_RT : array
            A 2d mesh of the Cartesian y coordinates.
        Z_RT : array
            A 3d array of yhe vaules for each point in x and y for each time.
        time : array
            The time array.
        constZ : float
            The constant z value (i.e. not the index)
        """
        #}}}

        self._X_RT   = X_RT
        self._Y_RT   = Y_RT
        self._Z_RT   = Z_RT
        self._time   = time
        self._constZ = constZ
    #}}}

    #{{{plotAndSavePerpPlane
    def plotAndSavePerpPlane(self):
        #{{{docstring
        """
        Performs the actual plotting of the perpendicular plane
        """
        #}}}

        # Initial plot (needed if we would like to save the plot)
        self._updatePerpAxInTime(0)

        # Call the save and show routine
        self.plotSaveAndShow(self._fig,\
                             self._updatePerpAxInTime,\
                             len(self._time))
    #}}}

    #{{{_updatePerpAxInTime
    def _updatePerpAxInTime(self, tInd):
        #{{{docstring
        """
        Updates the perpendicular axis.

        * Clears the axis
        * Plot the contourf
        * Set the labels
        * Updates the text
        * Updates the colorbar

        Parameters
        ----------
        tInd : int
            The current time index.
        """
        #}}}

        # Clear previous axis
        self._perpAx.cla()

        if self.__iterableLevels:
            self._cfKwargs.update({"vmax"   : self._vmax  [tInd],\
                                   "vmin"   : self._vmin  [tInd],\
                                   "levels" : self._levels[tInd],\
                                  })

        # Plot the perpoidal plane
        perpPlane = self._perpAx.\
            contourf(self._X_RT, self._Y_RT, self._Z_RT[tInd, :, :],\
                     **self._cfKwargs)

        # Set rasterization order
        self._perpAx.set_rasterization_zorder(self._axRasterization)
        # Draw the grids
        self._perpAx.grid(b=True)
        # Set x and y labels
        self._perpAx.\
            set_xlabel(self._ph.zTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)
        self._perpAx.\
            set_ylabel(self._ph.zTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)

        # Update the text
        self._updatePlotTxT(tInd)

        # Update the colorbar
        self._updateColorbar(self._fig, perpPlane, self._cBarAx)

        # Set equal axis
        self._perpAx.axis("equal")
    #}}}

    #{{{updatePlotTxt
    def updatePlotTxt(self, tInd):
        #{{{docstring
        """
        Updates the perpPlane plot by updating the axis, the title and
        formatting the axes

        Parameters
        ----------
        tInd : int
            The index to plot for
        perpAx : Axis
            The axis to update

        See the docstring of plotPerpPlane for details.
        """
        #}}}

        # Set title
        self._ph.zTxtDict["value"] = plotNumberFormatter(self._constVal, None)
        perpTitle = self._ph.zTxtDict["constZTxt"].format(self._ph.zTxtDict)
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._t[tInd], None, precision=4)
        timeTitle = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._perpAx.text(1.5, 1.05,\
                          "{}$,$ {}".format(perpTitle, timeTitle),\
                          transform = self._perpAx.transAxes,\
                          **self._txtKwargs)

        # Format axes
        self._perpAx.get_xaxis().\
            set_major_formatter(FuncFormatter(plotNumberFormatter))
        self._perpAx.get_yaxis().\
            set_major_formatter(FuncFormatter(plotNumberFormatter))
    #}}}
#}}}



# FIXME: Fix this!
#{{{updateParAx
def updateParAx(parAx, ax, X_RZ, Y_RZ, Z_RZ, Z_RZ_PPi,\
                ph, cfKwargs={}, labelSize=35):
    #{{{docstring
    """
    Updates the paroidal axis by cleaning it, and reset the plot

    This also sets the rasterization, the grid and the labels.

    Parameters
    ----------
    parAx : Axis
        The axis to update
    X_RZ : array
        A 2d mesh of the Cartesian x coordinates
    Y_RZ : array
        A 2d mesh of the Cartesian y coordinates
    Z_RZ : array
        The vaules for each point in x and y
    Z_RZ_PPi : array
        The vaules for each point in x and y, but pi away from the x and
        y of Z_RZ
    ph : PlotHelper
        The plot helper
    cfKwargs : dict
        Extra keywords for the contourf plot
    labelSize : int
        Size of the labelfont

    Returns
    -------
    parPlane : contour-like
        The object which can be used to find a colorbar
    """
    #}}}

    # Clear previous axis
    parAx.cla()

    # Plot the paroidal plane
    parPlane = ax.\
        contourf(X_RZ, Y_RZ, Z_RZ, **cfKwargs)

    # Also plot the negative plane
    ax.contourf(-X_RZ, Y_RZ, Z_RZ_PPi, **cfKwargs)

    # Set rasterization order
    parAx.set_rasterization_zorder(-10)
    # Draw the grids
    parAx.grid(b=True)
    # Set x and y labels
    parAx.set_xlabel(ph.zTxtDict["rhoTxtLabel"], fontsize = labelSize)
    parAx.set_ylabel(ph.zTxtDict["zTxtLabel"]  , fontsize = labelSize)

    return parPlane
#}}}

#{{{updatePolAx
def updatePolAx(polAx, ax, X_ZT, Y_ZT, Z_ZT, ph, cfKwargs={}, labelSize=35):
    #{{{docstring
    """
    Updates the poloidal axis by cleaning it, and reset the plot

    This also sets the rasterization, the grid and the labels.

    Parameters
    ----------
    polAx : Axis
        The axis to update
    X_ZT : array
        A 2d mesh of the Cartesian x coordinates
    Y_ZT : array
        A 2d mesh of the Cartesian y coordinates
    Z_ZT : array
        The vaules for each point in x and y
    ph : PlotHelper
        The plot helper
    cfKwargs : dict
        Extra keywords for the contourf plot
    labelSize : int
        Size of the labelfont

    Returns
    -------
    polPlane : contour-like
        The object which can be used to find a colorbar
    """
    #}}}

    # Clear previous axis
    polAx.cla()

    # Plot the poloidal plane
    polPlane = ax.\
        contourf(X_ZT, Y_ZT, Z_ZT.transpose(), **cfKwargs)

    # Set rasterization order
    polAx.set_rasterization_zorder(-10)
    # Draw the grids
    polAx.grid(b=True)
    # Set x and y labels
    polAx.set_xlabel(r"$\theta$", fontsize = labelSize)
    polAx.set_ylabel(ph.zTxtDict["zTxtLabel"], fontsize = labelSize)

    return polPlane
#}}}
