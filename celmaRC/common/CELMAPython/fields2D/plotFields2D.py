#!/usr/bin/env python

"""
Contains functions for plotting the 2D fields
"""

import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..superClasses import Plot2DSuperClass
from ..plotHelpers import plotNumberFormatter
import os


# Reset the size
plt.rc("xtick",  labelsize = 35)
plt.rc("ytick",  labelsize = 35)

#{{{Plot2DPerp
class Plot2DPerp(Plot2DSuperClass):
    """
    Class for 2D perpendicular plotting.

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
        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        self._fig = plt.figure(figsize = pltSize)
        self._fig.subplots_adjust(left=0.0)
        self._perpAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._perpAx).\
                 append_axes('right', '5%', '5%')
        self._perpAx.grid(True)
    #}}}

    #{{{setData
    def setData(self, X_RT, Y_RT, Z_RT, time, constZ, varName, savePath):
        #{{{docstring
        """
        Sets the data to be plotted (including the fileName and plot name)

        Parameters
        ----------
        X_RT : array
            A 2d mesh of the Cartesian x coordinates.
        Y_RT : array
            A 2d mesh of the Cartesian y coordinates.
        Z_RT : array
            A 3d array of the vaules for each point in x and y for each time.
        time : array
            The time array.
        constZ : float
            The constant z value (i.e. not the index).
        varName : str
            The name of the variable given in Z_RT.
        savePath : str
            Destination to save the plot in.
        """
        #}}}

        self._X_RT     = X_RT
        self._Y_RT     = Y_RT
        self._Z_RT     = Z_RT
        self._time     = time
        self._constZ   = constZ
        self._varName  = varName
        self._fileName =\
            os.path.join(savePath, "{}-{}-{}".format(varName, "perp", "2D"))

        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)
        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self._uc.conversionDict[self._varName])
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
        self.plotSaveShow(self._fig,\
                          self._fileName,\
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

        if self._iterableLevels:
            self._cfKwargs.update({"vmax"   : self._vmax  [tInd],\
                                   "vmin"   : self._vmin  [tInd],\
                                   "levels" : self._levels[tInd],\
                                  })

        # Plot the perpendicular plane
        perpPlane = self._perpAx.\
            contourf(self._X_RT, self._Y_RT, self._Z_RT[tInd, :, :],\
                     **self._cfKwargs)

        # Set rasterization order
        self._perpAx.set_rasterization_zorder(self._axRasterization)
        # Draw the grids
        self._perpAx.grid(b=True)
        # Set x and y labels
        self._perpAx.\
            set_xlabel(self._ph.rhoTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)
        self._perpAx.\
            set_ylabel(self._ph.rhoTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)

        # Update the text
        self._updatePlotTxt(tInd)

        # Update the colorbar
        self._updateColorbar(self._fig, perpPlane, self._cBarAx, tInd)

        # Set title
        self._cBar.set_label(self._varLabel, fontsize = self._labelSize)

        # Set equal axis
        self._perpAx.axis("equal")
    #}}}

    #{{{_updatePlotTxt
    def _updatePlotTxt(self, tInd):
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
        self._ph.zTxtDict["value"] = plotNumberFormatter(self._constZ, None)
        perpTitle = self._ph.zTxtDict["constZTxt"].format(self._ph.zTxtDict)
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[tInd], None, precision=4)
        timeTitle = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._perpAx.set_title("{}$,$ {}\n".format(perpTitle, timeTitle),\
                               fontsize = self._labelSize)

        # Format axes
        self._ph.makePlotPretty(self._perpAx,\
                                xprune   = "both",\
                                yprune   = "both",\
                                legend   = False,\
                                rotation = 20,\
                                )
    #}}}
#}}}

#{{{Plot2DPar
class Plot2DPar(Plot2DSuperClass):
    """
    Class for 2D parallel plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (10,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot2DPar

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
        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        self._fig = plt.figure(figsize = pltSize)
        self._fig.subplots_adjust(right=0.8)
        self._parAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._parAx).\
                 append_axes('right', '5%', '5%')
        self._parAx.grid(True)
    #}}}

    #{{{setData
    def setData(self, X_RZ, Y_RZ, Z_RZ, Z_RZ_PPi,\
                time, constTheta, varName, savePath):
        #{{{docstring
        """
        Sets the data to be plotted (including the fileName and plot name)

        Parameters
        ----------
        X_RZ : array
            A 2d mesh of the Cartesian x coordinates.
        Y_RZ : array
            A 2d mesh of the Cartesian y coordinates.
        Z_RZ : array
            A 3d array of the vaules for each point in x and y for each time.
        Z_RZ_PPi : array
            A 3d array of the vaules for each point in -x and y for each time.
        time : array
            The time array.
        constTheta : float
            The constant z value (i.e. not the index).
        varName : str
            The name of the variable given in Z_RZ.
        savePath : str
            Destination to save the plot in.
        """
        #}}}

        self._X_RZ       = X_RZ
        self._Y_RZ       = Y_RZ
        self._Z_RZ       = Z_RZ
        self._Z_RZ_PPi   = Z_RZ_PPi
        self._time       = time
        self._constTheta = int(constTheta)
        self._varName    = varName
        self._fileName   =\
            os.path.join(savePath, "{}-{}-{}".format(varName, "par", "2D"))

        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)
        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self._uc.conversionDict[self._varName])
    #}}}

    #{{{plotAndSaveParPlane
    def plotAndSaveParPlane(self):
        #{{{docstring
        """
        Performs the actual plotting of the parallel plane
        """
        #}}}

        # Initial plot (needed if we would like to save the plot)
        self._updateParAxInTime(0)

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updateParAxInTime,\
                          len(self._time))
    #}}}

    #{{{_updateParAxInTime
    def _updateParAxInTime(self, tInd):
        #{{{docstring
        """
        Updates the parallel axis.

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
        self._parAx.cla()

        if self._iterableLevels:
            self._cfKwargs.update({"vmax"   : self._vmax  [tInd],\
                                   "vmin"   : self._vmin  [tInd],\
                                   "levels" : self._levels[tInd],\
                                  })

        # Plot the parallel plane
        parPlane = self._parAx.\
            contourf(self._X_RZ, self._Y_RZ, self._Z_RZ[tInd, :, :],\
                     **self._cfKwargs)
        # Plot the negative parallel plane
        parPlane =\
        self._parAx.\
            contourf(-self._X_RZ, self._Y_RZ, self._Z_RZ_PPi[tInd, :, :],\
                     **self._cfKwargs)

        # Set rasterization order
        self._parAx.set_rasterization_zorder(self._axRasterization)
        # Draw the grids
        self._parAx.grid(b=True)
        # Set x and y labels
        self._parAx.\
            set_xlabel(self._ph.rhoTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)
        self._parAx.\
            set_ylabel(self._ph.zTxtDict["zTxtLabel"],\
                       fontsize = self._labelSize)

        # Update the text
        self._updatePlotTxt(tInd)

        # Update the colorbar
        self._updateColorbar(self._fig, parPlane, self._cBarAx, tInd)

        # Set title
        self._cBar.set_label(self._varLabel, fontsize = self._labelSize)
    #}}}

    #{{{_updatePlotTxt
    def _updatePlotTxt(self, tInd):
        #{{{docstring
        """
        Updates the parPlane plot by updating the axis, the title and
        formatting the axes

        Parameters
        ----------
        tInd : int
            The index to plot for
        parAx : Axis
            The axis to update

        See the docstring of plotParPlane for details.
        """
        #}}}

        # Set title
        self._ph.thetaTxtDict["value"] = plotNumberFormatter(self._constTheta, None)
        parTitle = self._ph.thetaTxtDict["constThetaTxt"].format(self._ph.thetaTxtDict)
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[tInd], None, precision=4)
        timeTitle = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._parAx.set_title("{}$,$ {}\n".format(parTitle, timeTitle),\
                               fontsize = self._labelSize)

        # Format axes
        self._ph.makePlotPretty(self._parAx,\
                                xprune   = "both",\
                                yprune   = "both",\
                                legend   = False,\
                                rotation = 20,\
                                )
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
