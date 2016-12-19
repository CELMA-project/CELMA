#!/usr/bin/env python

"""
Contains super class for plotting the 2D fields
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
                 convertToPhysical,\
                 show = False,\
                 save = True,\
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
        * Sets the var label template

        Parameters
        ----------
        path : str
            Destination of save
        fluct: bool
            Whether or not the fluctuations are being plotted.
        convertToPhysical : bool
            Whether or not to convert to physical units.
        show : bool
            Whether or not the plot is to be displayed.
        save : bool
            Whether or not to save the plot.
        savePath : str
            Path (excluding extension) to use if the plot or animation is
            saved.
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
        self.setCollectPath(path, convertToPhysical, uc)

        # Set animation and text options
        self.setAnimationOptions()

        # Set extra contourf keyword arguments
        self._cfKwargs = {}
        self.setContourfArguments()
        if self._fluct:
            cmap = divCMap
        else:
            cmap = seqCMap

        # Update extra contourf keyword arguments with cmap and zorder
        self._cfKwargs.update({"cmap" : cmap, "zorder" : -20})

        # Magic numbers
        self._labelSize = 35
        self._axRasterization = -10

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

    #{{{setCollectPath
    def setCollectPath(self, path, convertToPhysical=True, uc=None):
        #{{{docstring
        """
        Sets the collect path, makes the units converter and plot helper

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
        self.collectPath = path

        if uc is None:
            self._uc = UnitsConverter(path, convertToPhysical)

        # Make the plot helper
        self._ph = PlotHelper(convertToPhysical)
        self._ph.makeDimensionStringsDicts(self._uc)
        self._convertToPhysical = self._uc.convertToPhysical
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

            if self._save:
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
            if self._save:
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

        if self._show:
            fig.show()

        plt.close(fig)
    #}}}
#}}}
