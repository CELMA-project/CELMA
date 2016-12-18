#!/usr/bin/env python

"""
Contains functions for plotting the 2D fields
"""

from ..superClasses import Plot2DSuperClass
from ..plotHelpers import plotNumberFormatter
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pylab as plt
import numpy as np
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
    def __init__(self, *args, pltSize = (20,15), **kwargs):
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

        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        self._fig = plt.figure(figsize = pltSize)
        self._fig.subplots_adjust(left=0.0)
        self._perpAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._perpAx).\
                 append_axes('right', '5%', '5%')
        self._perpAx.grid(True)

        # Set the axis title
        self._axTitle = "{}$,$ {}\n"
    #}}}

    #{{{setPerpData
    def setPerpData(self, X_RT, Y_RT, Z_RT, time, constZ, varName, savePath):
        #{{{docstring
        """
        Sets the perpendicular data and label to be plotted

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
        self._savePath = savePath

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

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._varName, "perp", "2D"))

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
        self._updatePerpPlotTxt(tInd)

        # Update the colorbar
        self._updateColorbar(self._fig, perpPlane, self._cBarAx, tInd)

        # Set title
        self._cBar.set_label(self._varLabel, fontsize = self._labelSize)

        # Set equal axis
        self._perpAx.axis("equal")
    #}}}

    #{{{_updatePerpPlotTxt
    def _updatePerpPlotTxt(self, tInd):
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
        self._perpAx.set_title(self._axTitle.format(perpTitle, timeTitle),\
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
    def __init__(self, *args, pltSize = (20,15), **kwargs):
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
        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        self._fig = plt.figure(figsize = pltSize)
        self._fig.subplots_adjust(right=0.8)
        self._parAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._parAx).\
                 append_axes('right', '5%', '5%')
        self._parAx.grid(True)

        # Set the axis title
        self._axTitle = "{}$,$ {}\n"
    #}}}

    #{{{setParData
    def setParData(self, X_RZ, Y_RZ, Z_RZ, Z_RZ_PPi,\
                time, constTheta, varName, savePath):
        #{{{docstring
        """
        Sets the parallel data and label to be plotted

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
        self._savePath   = savePath

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

        # Set the file name
        self._fileName   =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._varName, "par", "2D"))

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
        self._updateParPlotTxt(tInd)

        # Update the colorbar
        self._updateColorbar(self._fig, parPlane, self._cBarAx, tInd)

        # Set title
        self._cBar.set_label(self._varLabel, fontsize = self._labelSize)
    #}}}

    #{{{_updateParPlotTxt
    def _updateParPlotTxt(self, tInd):
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
        self._parAx.set_title(self._axTitle.format(parTitle, timeTitle),\
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

#{{{Plot2DPol
class Plot2DPol(Plot2DSuperClass):
    """
    Class for 2D poloidal plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (20,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot2DPol

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
        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        self._fig = plt.figure(figsize = pltSize)
        self._fig.subplots_adjust(right=0.8)
        self._polAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._polAx).\
                 append_axes('right', '5%', '5%')
        self._polAx.grid(True)

        # Set the axis title
        self._axTitle = "{}$,$ {}\n"
    #}}}

    #{{{setData
    def setData(self, X_ZT, Y_ZT, Z_ZT, time, constRho, varName, savePath):
        #{{{docstring
        """
        Sets the poloidal data and label to be plotted

        Parameters
        ----------
        X_ZT : array
            A 2d mesh of the Cartesian x coordinates.
        Y_ZT : array
            A 2d mesh of the Cartesian y coordinates.
        Z_ZT : array
            A 3d array of the vaules for each point in x and y for each time.
        time : array
            The time array.
        constRho : float
            The constant z value (i.e. not the index).
        varName : str
            The name of the variable given in Z_ZT.
        savePath : str
            Destination to save the plot in.
        """
        #}}}

        self._X_ZT     = X_ZT
        self._Y_ZT     = Y_ZT
        self._Z_ZT     = Z_ZT
        self._time     = time
        self._constRho = constRho
        self._varName  = varName

        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)
        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self._uc.conversionDict[self._varName])
    #}}}

    #{{{plotAndSavePolPlane
    def plotAndSavePolPlane(self):
        #{{{docstring
        """
        Performs the actual plotting of the poloidal plane
        """
        #}}}

        # Set the file name
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._varName, "pol", "2D"))

        # Initial plot (needed if we would like to save the plot)
        self._updatePolAxInTime(0)

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updatePolAxInTime,\
                          len(self._time))
    #}}}

    #{{{_updatePolAxInTime
    def _updatePolAxInTime(self, tInd):
        #{{{docstring
        """
        Updates the poloidal axis.

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
        self._polAx.cla()

        if self._iterableLevels:
            self._cfKwargs.update({"vmax"   : self._vmax  [tInd],\
                                   "vmin"   : self._vmin  [tInd],\
                                   "levels" : self._levels[tInd],\
                                  })

        # Plot the poloidal plane
        polPlane = self._polAx.\
            contourf(self._X_ZT, self._Y_ZT, self._Z_ZT[tInd, :, :].transpose(),\
                     **self._cfKwargs)

        # Set rasterization order
        self._polAx.set_rasterization_zorder(self._axRasterization)
        # Draw the grids
        self._polAx.grid(b=True)
        # Set x and y labels
        self._polAx. set_xlabel(r"$\theta$", fontsize = self._labelSize)
        self._polAx.\
            set_ylabel(self._ph.zTxtDict["zTxtLabel"],\
                       fontsize = self._labelSize)

        # Update the text
        self._updatePolPlotTxt(tInd)

        # Tweak latex on x-axis
        self._polAx.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
        self._polAx.set_xticklabels(\
                (r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"))

        # Update the colorbar
        self._updateColorbar(self._fig, polPlane, self._cBarAx, tInd)

        # Set title
        self._cBar.set_label(self._varLabel, fontsize = self._labelSize)
    #}}}

    #{{{_updatePolPlotTxt
    def _updatePolPlotTxt(self, tInd):
        #{{{docstring
        """
        Updates the polPlane plot by updating the axis, the title and
        formatting the axes

        Parameters
        ----------
        tInd : int
            The index to plot for
        polAx : Axis
            The axis to update

        See the docstring of plotPolPlane for details.
        """
        #}}}

        # Set title
        self._ph.rhoTxtDict["value"] = plotNumberFormatter(self._constRho, None)
        polTitle = self._ph.rhoTxtDict["constRhoTxt"].format(self._ph.rhoTxtDict)
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[tInd], None, precision=4)
        timeTitle = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._polAx.set_title(self._axTitle.format(polTitle, timeTitle),\
                               fontsize = self._labelSize)

        # Format axes
        self._ph.makePlotPretty(self._polAx,\
                                xprune   = "both",\
                                yprune   = "both",\
                                legend   = False,\
                                )
    #}}}
#}}}

#{{{Plot2DPerpPar
class Plot2DPerpPar(Plot2DPerp, Plot2DPar):
    """
    Class for 2D perpendicular-parallel plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (35,15), **kwargs):
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

        # Re-open the figure
        plt.close("all")
        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        self._fig = plt.figure(figsize = pltSize)

        # Create figure and axes
        # NOTE: tight_layout=True gives wobbly plot as the precision of
        #       the colorbar changes during the animation
        gs           = GridSpec(1, 3, width_ratios=(20, 20, 1))
        self._perpAx = self._fig.add_subplot(gs[0])
        self._parAx  = self._fig.add_subplot(gs[1])
        self._cBarAx = self._fig.add_subplot(gs[2])
        self._fig.subplots_adjust(wspace=0.25)
        self._parAx.grid(True)
        self._perpAx.grid(True)

        # Set the axis title
        self._axTitle = "{}\n"
    #}}}

    #{{{plotAndSavePerpPlane
    def plotAndSavePerpPlane(self):
        #{{{docstring
        """
        Performs the actual plotting of the perpendicular and parallel
        """
        #}}}

        # Set the file name
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._varName, "perpPar", "2D"))

        # Initial plot (needed if we would like to save the plot)
        self._updatePerpAndPolAxInTime(0)

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updatePerpAndPolAxInTime,\
                          len(self._time))
    #}}}

    #{{{_updatePerpAndPolAxInTime
    def _updatePerpAndPolAxInTime(self, tInd):
        #{{{docstring
        """
        Updates the perpendicular and parallel axis and sets the figure title.

        Parameters
        ----------
        tInd : int
            The current time index.
        """
        #}}}

        self._updatePerpAxInTime(tInd)
        self._updateParAxInTime(tInd)

        timeTitle = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}\n\n\n".\
                           format(timeTitle), fontsize = self._labelSize)
    #}}}
#}}}
