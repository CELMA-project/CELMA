#!/usr/bin/env python

"""
Contains functions for plotting the 2D fields
"""

from ..plotHelpers import (PlotHelper,\
                           plotNumberFormatter,\
                           seqCMap,\
                           divCMap)
from ..unitsConverter import UnitsConverter
import numpy as np
import matplotlib.animation as animation
import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
import os

#{{{Plot2DSuperClass
class Plot2DSuperClass(object):
    """
    Super class for 2D plots.

    Handles common plot options and saving.
    """

    #{{{constructor
    def __init__(self,\
                 path,\
                 fluct,\
                 show = False,\
                 save = True,\
                 convertToPhysical = True,\
                 extension = None,\
                 uc = None):
        #{{{docstring
        """
        Constructor for Plot2DSuperClass.

        * Stores common plotting options:
            * Animation
            * Text
            * Extra contourf arguments
        * Makes the unitsConverter if not given as input
        * Makes the plot helper object

        Parameters
        ----------
        path : str
            Destination of save
        fluct: bool
            Whether or not the fluctuations are being plotted.
        show : bool
            Whether or not the plot is to be displayed.
        save : bool
            Whether or not to save the plot.
        savePath : str
            Path (excluding extension) to use if the plot or animation is
            saved.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        extension : [None|str]
            Overrides default extension if set. Excludes the "."
        uc : [None|UnitsConverter]
            The UnitsConverter will be set from the path if not given.
        """
        #}}}

        # Set memberdata
        self._fluct     = fluct
        self._show      = show
        self._save      = save
        self._extension = extension

        # Set the path, units converter and plot helper
        self.setPath(path, convertToPhysical, uc)

        # Set animation and text options
        self.setAnimationOptions()
        self.setTxTProperties()

        # Set extra contourf keyword arguments
        self.setContourfArguments()
        if self._fluct:
            cmap = divCMap
        else:
            cmap = seqCMap

        # Update extra contourf keyword arguments with cmap and zorder
        self._cfKwargs.update({"cmap" : cmap, "zorder" : -20})

        # Magic numbers
        self._labelSize=35
        self._axRasterization-10
    #}}}

    #{{{_updateColorbar
    def _updateColorbar(self, fig, plane, cBarAx):
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
        """
        #}}}

        # Set the colorbar
        # Clear the axis
        # http://stackoverflow.com/questions/39472017/how-to-animate-the-colorbar-in-matplotlib/39596853
        cBarAx.cla()
        if self._fluct:
            # Create the ticks (11 with 0 in the center)
            nTicks = 11
            ticks  = np.linspace(self._varMin, self._varMax, nTicks)
            # Enforce the center one to be 0 (without round off)
            ticks[int((nTicks - 1)/2)] = 0
        else:
            ticks = None
        fig.colorbar(plane,\
                     cax    = cBarAx,\
                     ticks  = ticks,\
                     format = FuncFormatter(plotNumberFormatter))
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

    #{{{setPath
    def setPath(self, path, convertToPhysical=True, uc=None):
        #{{{docstring
        """
        Sets the path, makes the units converter and plot helper

        Parameters
        ----------
        path : str
            Destination of save
        convertToPhysical : bool
            Whether or not to convert to physical units.
        uc : [None|UnitsConverter]
            The UnitsConverter will be set from the path if not given.
        """
        #}}}
        self._path = path

        if uc is None:
            self._uc = UnitsConverter(path, convertToPhysical)

        # Make the plot helper
        self._ph = PlotHelper(convertToPhysical)
        self._convertToPhysical = self._ph.convertToPhysical
    #}}}

    #{{{setTxTProperties
    def setTxTProperties(self, ha="center", va="center", fontsize=35):
        #{{{docstring
        """
        Set the text properties for the titles

        Parameters
        ----------
        ha : str
            Horizontal alignment of the text
        va : str
            Vertical algnemnt of the text

        """
        #}}}
        self._txtKwargs = {"ha"       : ha      ,\
                           "va"       : va      ,\
                           "fontsize" : fontsize,\
                          }
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
        if vmax == None:
            self._iterableLevels = False
            if vmin is not None and levels is not None:
                success = False
        elif vmin == None:
            self._iterableLevels = False
            if vmax is not None and levels is not None:
                success = False
        elif levels == None:
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

            self._cfKwargs = {"vmax"   : self._varMax[0],\
                              "vmin"   : self._varMin[0],\
                              "levels" : self._levels[0],\
                             }
        else:
            self._cfKwargs = {"vmax"   : None,\
                              "vmin"   : None,\
                              "levels" : None,\
                             }
    #}}}

    #{{{plotSaveAndShow
    def plotSaveAndShow(self, fig, fileName, func, fargs, frames):
        #{{{docstring
        """
        Saves and closes the plot.

        Parameters
        ----------
        fig : Figure
            Figure to save.
        fileName : str
            Name of the file, excluding the path and the extension
        func : function
            The function to use for generating the animation.
        fargs : tuple
            Tuple of the arguments to use in the animation.
        frames : int
            Number of frames.
            If this is less than one, a normal plot will be made.
        """
        #}}}

        # Merge fileName and path
        fileName = os.path.join(self._path, fileName)
        if frames > 1:
            # Animate
            anim = animation.FuncAnimation(fig            ,\
                                           func           ,\
                                           fargs  = fargs ,\
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