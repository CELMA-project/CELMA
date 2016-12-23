#!/usr/bin/env python

"""
Contains super classes for animations plots
"""

from ..plotHelpers import (PlotHelper,\
                           plotNumberFormatter,\
                           seqCMap,\
                           seqCMap3,\
                           divCMap)
from ..unitsConverter import UnitsConverter
from ..plotHelpers import getMaxMinAnimation
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
import numpy as np
import matplotlib.animation as animation
import matplotlib.pylab as plt
import os

#{{{PlotAnimSuperClass
class PlotAnimSuperClass(PlotSuperClass):
    """
    Super class for animation plots.

    Handles common plot options and saving.
    """

    #{{{constructor
    def __init__(self   ,\
                 *args  ,\
                 *kwargs):
        #{{{docstring
        """
        Constructor for PlotAnimSuperClass.

        * Calls the parent constructor
        * Sets the animation options

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

        # Set animation and text options
        self.setAnimationOptions()

        # Magic numbers
        self._labelSize = 35
        self._axRasterization = -10
    #}}}

    #{{{setAnimationOptions
    def setAnimationOptions(self, bitrate = -1, fps = 10, codec = "h264"):
        #{{{docstring
        """
        Reset the save destination.

        Parameters
        ----------
        bitrate : int
            Sets the quality of the animation.
            -1 gives automatic quality.
        fps : int
            Frames per second.
            Sets the speed of the plot as indicated in
            http://stackoverflow.com/questions/22010586/matplotlib-animation-duration
        codec : str
            Codec to use.
            Default is h264.
            For installation, see
            https://github.com/loeiten/usingLinux/blob/master/installationProcedures/ffmpeg.md
        """
        #}}}
        self._bitrate = bitrate
        self._fps     = fps
        self._codec   = codec
    #}}}

    #{{{plotSaveShow
    def plotSaveShow(self, fig, fileName, func, frames):
        #{{{docstring
        """
        Saves, show and closes the plot.

        Parameters
        ----------
        fig : Figure
            Figure to save.
        fileName : str
            Name of the file, including the path and the excluding the
            extension
        func : function
            The function to use for generating the animation.
        frames : int
            Number of frames.
            If this is less than one, a normal plot will be made.
        """
        #}}}

        if frames > 1:
            # Animate
            anim = animation.FuncAnimation(fig            ,\
                                           func           ,\
                                           frames = frames,\
                                           blit   = False ,\
                                           )

            if self._savePlot:
                FFMpegWriter = animation.writers['ffmpeg']
                writer = FFMpegWriter(bitrate = self._bitrate,\
                                      fps     = self._fps    ,\
                                      codec   = self._codec)

                if self._extension is None:
                    self._extension = "mp4"

                # Save the animation
                fileName = "{}.{}".format(fileName, self._extension)
                anim.save(fileName, writer = writer)
                print("Saved to {}".format(fileName))
        else:
            if self._savePlot:
                if self._extension is None:
                    self._extension = "png"

                # Save the figure
                fileName = "{}.{}".format(fileName, self._extension)
                fig.savefig(fileName,\
                            transparent = True  ,\
                            bbox_inches = "tight",\
                            pad_inches  = 0      ,\
                            )
                print("Saved to {}".format(fileName))

        if self._showPlot:
            fig.show()

        plt.close(fig)
    #}}}
#}}}

#{{{PlotAnim1DSuperClass
class PlotAnim1DSuperClass(PlotAnimSuperClass):
    """
    Super class for 1D animation plots.

    Handles common plot options and saving.
    """

    #{{{constructor
    def __init__(self,\
                 *args,\
                 **kwargs):
        #{{{docstring
        """
        Constructor for PlotAnim1DSuperClass.

        * Calls the parent constructor
        * Sets the variable legend template

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

        # Set the colormap
        self._cmap = seqCMap3

        # Magic numbers
        self._cols     = 2
        self._marker   = "o"
        self._ddtColor = "k"

        # Set var label template
        if self._convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
        else:
            unitsOrNormalization = "${normalization}$"

        self._varLegendTemplate = r"${{}}${}".format(unitsOrNormalization)
    #}}}

    #{{{_createFiguresAndAxes
    def _createFiguresAndAxes(self):
        """
        Creates the figures and axes and organizes the plot.
        """

        # Calculate the number of rows
        rows = int(np.ceil(len(self._vars)/self._cols))

        # Create the figure
        self._fig = plt.figure(figsize=self._pltSize)

        # Create the axes
        gs = GridSpec(rows, self._cols)
        self._axes = []
        # Make first ax
        firstAx = self._fig.add_subplot(gs[0])
        firstAx.grid(True)
        self._axes.append(firstAx)
        # Make the rest of the axes
        for nr in range(1, len(self._vars)):
            ax = self._fig.add_subplot(gs[nr])
            ax.grid(True)
            self._axes.append(ax)

        # Cast to tuple
        self._axes = tuple(self._axes)

        # Set the x-labels
        totalAx = len(tuple(self._vars.keys()))
        for nr in range(totalAx-self._cols):
            self._axes[nr].tick_params(labelbottom="off")
        for nr in range(totalAx-self._cols, totalAx):
            self._axes[nr].set_xlabel(self._xlabel)

        # Adjust the subplots
        self._fig.subplots_adjust(hspace=0.1, wspace=0.35)
    #}}}

    #{{{_initialPlot
    def _initialPlot(self):
        #{{{docstring
        """
        Initial plot.

        The initial plot:
            * Fetches the labels
            * Finds the max/min of all lines
            * Plots normal lines and makes the axes pretty
            * Plots ddt (if present)
            * Populates self._lines and self._ddtLines with the lines
        """
        #}}}

        # Placeholder for the lines
        self._lines = []

        # Plot and obtain the level
        for ax, key, color in zip(self._axes, self._plotOrder, self._colors):
            # Obtain the legend
            pltVarName = self._ph.getVarPltName(key)
            if self._convertToPhysical:
                legend = self._varLegendTemplate.\
                    format(pltVarName, **self._uc.conversionDict[key])
            else:
                legend = "${}$".format(pltVarName)
            # Do the plot
            self._lines.append(\
                        ax.plot(self._X,\
                                self._vars[key][0,:],\
                                marker          = self._marker,\
                                color           = color       ,\
                                markeredgecolor = color       ,\
                                markerfacecolor = color       ,\
                                label           = legend      ,\
                                )[0]\
                              )

            # FIXME: varyMaxMin currently not implemented
            vMax, vMin = getMaxMinAnimation((self._vars[key],), False, False)
            ax.set_ylim(vMin[0], vMax[0])

            # Make the ax pretty
            self._ph.makePlotPretty(ax,\
                                    yprune="both",\
                                    loc="upper right",\
                                    rotation = 45,\
                                    ybins = 6
                                    )

        # If ddt is present
        if self._ddtPresent:
            self._ddtLines = []
            for key, color in zip(self._plotOrder, self._colors):
                # Redo the plots
                self._ddtLines.append(\
                        self._axes[-1].plot(self._X                       ,\
                                            self._vars[key][0,:]          ,\
                                            marker          = self._marker,\
                                            color           = color       ,\
                                            markeredgecolor = color       ,\
                                            markerfacecolor = color       ,\
                                            )[0]\
                                  )
            pltVarName = self._ph.getVarPltName(self._ddtVar)
            if self._convertToPhysical:
                legend = self._varLegendTemplate.\
                    format(pltVarName, **self._uc.conversionDict[key])
            else:
                legend = "${}$".format(pltVarName)

            # Plot the ddt
            self._ddtLines.append(\
                    self._axes[-1].plot(self._X                         ,\
                                        self._vars[self._ddtVar][0,:]   ,\
                                        marker          = self._marker  ,\
                                        color           = self._ddtColor,\
                                        markeredgecolor = self._ddtColor,\
                                        markerfacecolor = self._ddtColor,\
                                        label           = legend        ,\
                                        )[0]\
                                )

            # Make an array tuple in order to check for max/min
            arrayTuple = tuple(self._vars[key] for key in self._plotOrder)
            # FIXME: varyMaxMin currently not implemented
            vMax, vMin = getMaxMinAnimation(arrayTuple, False, False)
            self._axes[-1].set_ylim(vMin[0], vMax[0])

            # Make the ax pretty
            self._ph.makePlotPretty(self._axes[-1],\
                                    yprune="both",\
                                    loc="upper right",\
                                    rotation = 45,\
                                    ybins = 6)

        # Cast to tuples
        self._lines = tuple(self._lines)
        if self._ddtPresent:
            self._ddtLines = tuple(self._ddtLines)
    #}}}

    #{{{_setColors
    def _setColors(self):
        """
        Sets the colors to be used in the plotting
        """
        colorSpace = np.arange(len(self._vars))
        self._colors = self._cmap(np.linspace(0, 1, len(colorSpace)))
    #}}}
#}}}

#{{{PlotAnim2DSuperClass
class PlotAnim2DSuperClass(PlotAnimSuperClass):
    """
    Super class for 2D animation plots.

    Handles common plot options and saving.
    """

    #{{{constructor
    def __init__(self,\
                 *args,\
                 fluct = None,\
                 **kwargs):
        #{{{docstring
        """
        Constructor for PlotAnim2DSuperClass.

        * Calls the parent constructor

        * Stores common plotting options:
            * Text
            * Extra contourf arguments
        * Sets the var label template

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details
        fluct: bool
            Whether or not the fluctuations are being plotted.
        **kwargs : keyword arguments
            See parent constructor for details
        """
        #}}}

        # Guard
        if fluct is None:
            raise ValueError("'fluct' must be bool")

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set memberdata
        self._fluct = fluct

        # Set extra contourf keyword arguments
        self._cfKwargs = {}
        self.setContourfArguments()

        # Update extra contourf keyword arguments with cmap and zorder
        if self._fluct:
            cmap = divCMap
        else:
            cmap = seqCMap
        self._cfKwargs.update({"cmap" : cmap, "zorder" : -20})

        # Set var label template
        if self._convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
        else:
            unitsOrNormalization = "${normalization}$"
        if not(self._fluct):
            self._varLabelTemplate = r"${{}}${}".format(unitsOrNormalization)
        else:
            self._varLabelTemplate = r"$\widetilde{{{{{{}}}}}}${}".\
                    format(unitsOrNormalization)
    #}}}

    #{{{_updateColorbar
    def _updateColorbar(self, fig, plane, cBarAx, tInd):
        #{{{docstring
        """
        Updates the colorbar.

        Parameters
        ----------
        fig : Figure
            The figure where the colorbar belongs.
        plane : contour-like
            The plane where the colorbar reads the data from.
        cBarAx : Axis
            The axis where the colorbar belongs.
        tInd : int
            The current time index.
        """
        #}}}

        # Set the colorbar
        # Clear the axis
        # http://stackoverflow.com/questions/39472017/how-to-animate-the-colorbar-in-matplotlib/39596853
        cBarAx.cla()
        if self._fluct and self._iterableLevels:
            # Create the ticks (11 with 0 in the center)
            nTicks = 11
            ticks  = np.linspace(self._vmin[tInd], self._vmax[tInd], nTicks)
            # Enforce the center one to be 0 (without round off)
            ticks[int((nTicks - 1)/2)] = 0
        else:
            ticks = None

        self._cBar = fig.colorbar(plane,\
                        cax    = cBarAx,\
                        ticks  = ticks,\
                        format = FuncFormatter(plotNumberFormatter))
    #}}}

    #{{{setContourfArguments
    def setContourfArguments(self, vmax=None, vmin=None, levels=None):
        #{{{docstring
        """
        Set extra contourf keyword arguments.

        Parameters
        ---------
        vmax : [None|tuple]
            Max value to give a color.
            One for each frame if not None.
            If None: vmin and levels must also be None
            See getMaxMinAnimation in plot2DHelpers for a function which
            automatically does this.
        vmin : [None|tuple]
            Min value to give a color.
            One for each frame if not None.
            If None: vmax and levels must also be None
            See getMaxMinAnimation in plot2DHelpers for a function which
            automatically does this.
        levels : [None|tuple]
            The levels to  use.
            One for each frame if not None.
            If None: vmax and vmin must also be None
            See getLevelsAnimation in plot2DHelpers for a function which
            automatically does this.
        """
        #}}}

        # Guard
        success = True
        self._iterableLevels = True
        if vmax is None:
            self._iterableLevels = False
            if vmin is not None and levels is not None:
                success = False
        elif vmin is None:
            self._iterableLevels = False
            if vmax is not None and levels is not None:
                success = False
        elif levels is None:
            self._iterableLevels = False
            if vmax is not None and vmin is not None:
                success = False
        if not(success):
            message = "Either all or none of vmax, vmin and levels must be None"
            raise ValueError(message)

        if self._iterableLevels:
            self._vmax   = vmax
            self._vmin   = vmin
            self._levels = levels

            self._cfKwargs.update({"vmax"   : self._vmax[0]  ,\
                                   "vmin"   : self._vmin[0]  ,\
                                   "levels" : self._levels[0],\
                                   })
        else:
            self._cfKwargs.update({"vmax"   : None,\
                                   "vmin"   : None,\
                                   "levels" : None,\
                                  })
    #}}}

    #{{{_setFileName
    def _setFileName(self, plotTypeName):
        #{{{docstring
        """
        Sets self._fileName

        Parameters
        ----------
        plotTypeName : str
            The plot type name
        """
        #}}}
        if self._fluct:

            fileName = "{}-{}-{}-fluct"\
                    .format(self._varName, plotTypeName, "2D")
        else:
            fileName = "{}-{}-{}".format(self._varName, plotTypeName, "2D")

        self._fileName = os.path.join(self._savePath, fileName)
    #}}}
#}}}
