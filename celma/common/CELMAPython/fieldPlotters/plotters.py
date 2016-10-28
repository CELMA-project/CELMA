#!/usr/bin/env python

"""
Contains classes for plotting the fields
"""

from ..plotHelpers import (PlotHelper,\
                           plotNumberFormatter,\
                           seqCMap,\
                           divCMap,\
                           findLargestRadialGrad)
from ..statsAndSignals import polAvg
from .cylinderMesh import CylinderMesh
from matplotlib import get_backend
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from boutdata import collect
import numpy as np
import os
import warnings

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called "folder" in
# __call_post_processing_function)

#{{{class Plot
class Plot(object):
    """
    Parent class for plotting the results of the CELMA code.

    Handles:

    * Setting the collect options
    * Formation of ticks
    * Preparation of xlabel, ylabel and title
    """

    #{{{Constructor
    def __init__(self                             ,\
                 path                             ,\
                 xguards           = False        ,\
                 yguards           = False        ,\
                 xSlice            = slice(0,None),\
                 ySlice            = slice(0,None),\
                 zSlice            = slice(0,None),\
                 tSlice            = None         ,\
                 convertToPhysical = False        ,\
                 subPolAvg         = False        ,\
                 showPlot          = False        ,\
                 savePlot          = True         ,\
                 saveFolder        = None         ,\
                 extension         = "png"        ,\
                 writer            = "ffmpeg"     ,\
                ):
        #{{{docstring
        """
        The constructor:

        1. Sets the plot style
        2. Calcuates rho, theta and z
        3. Collects the time
        4. Collects normalizing parameters if set

        Parameters
        ----------
        path : str
            The path to collect from.
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        xSlice : slice
            How the data will be sliced in x.
        ySlice : slice
            How the data will be sliced in y.
        zSlice : slice
            How the data will be sliced in z.
        tSlice : slice
            How the data will be sliced in t.
        maxGradRhoFolder : [None | str]
            If this is set, the xSlice will be replaced by the index of
            the largest gradient in rho direction.
        convertToPhysical : bool
            If the physical or normalized units should be plotted.
        subPolAvg : vool
            Whether or not the poloidal average should be.
            subtracted from the data.
        showPlot : bool
            If the plot should be displayed.
        savePlot : bool
            If plot should be saved.
        saveFolder : str
            Name of the folder to save plots in.
        extension : str
            Extension of the plot (if the animation is not used).
        writer : str
            Writer to use if the plots are animated.
        """
        #}}}

        # Set the bitrates, fps and codec for ffmpeg (currently magic numbers)
        self._bitrate = -1
        self._fps     = 10
        self._codec   = "h264"

        # Set member data from input
        self._path       = path
        self._xguards    = xguards
        self._yguards    = yguards
        self._showPlot   = showPlot
        self._savePlot   = savePlot
        self._saveFolder = saveFolder
        self._subPolAvg  = subPolAvg
        self._extension  = extension
        self._writer     = writer

        if writer == "imagemagick":
            self._animExtension = ".gif"
        elif writer == "ffmpeg":
            self._animExtension = ".mp4"
        else:
            message = "Writer {} not implemented".format(writer)
            raise NotImplementedError(message)

        # Get proper indices
        self._xind = self._getIndices(xSlice, "x")
        self._yind = self._getIndices(ySlice, "y")
        self._zind = self._getIndices(zSlice, "z")
        self._tind = self._getIndices(tSlice, "t")

        if type(xSlice) == slice:
            if (xSlice.step is not None):
                message = ("{0}{1}WARNING: xSlice.step not implemented.\n"
                           "Setting to None{1}{0}".format("\n"*3, "!"*3))
                print(message)
                xSlice.step = None
        if type(ySlice) == slice:
            if (ySlice.step is not None):
                message = ("{0}{1}WARNING: ySlice.step not implemented.\n"
                           "Setting to None{1}{0}".format("\n"*3, "!"*3))
                print(message)
                ySlice.step = None
        if type(zSlice) == slice:
            if (zSlice.step is not None):
                message = ("{0}{1}WARNING: ySlice.step not implemented.\n"
                           "Setting to None{1}{0}".format("\n"*3, "!"*3))
                print(message)
                zSlice.step = None

        # Used if we are taking poloidal averages
        self._xSlice = xSlice
        self._ySlice = ySlice
        self._zSlice = zSlice
        # Used to chop the data
        self._tSlice = tSlice

        # Get the time
        t = collect("t_array", path=self._path, tind=self._tind, info=False)

        # Slice in t
        if self._tSlice is not None:
            if self._tSlice.step is not None:
                t = t[::self._tSlice.step]

        # Set frames
        self._frames = len(t)

        # Make the PlotHelper object
        # Public as used in the driver
        self.helper = PlotHelper(path                                 ,\
                                 t                 = t                ,\
                                 xguards           = xguards          ,\
                                 yguards           = yguards          ,\
                                 convertToPhysical = convertToPhysical,\
                                 )

        # Set colormap
        if self._subPolAvg:
            self._cmap = divCMap
        else:
            self._cmap = seqCMap
    #}}}

    #{{{ _getIndices
    def _getIndices(self, curSlice, dimension):
        """
        Return the slice such that it can be given as an input to "collect"

        Parameters
        ----------
        curSlice : [slice | int | None]
            Current slice to use
        dimension : ["x" | "y" | "z" | "t"]
            The dimension to slice in

        Returns
        -------
        curIndices : sequence (not str)
            A sequence of the start and the stop values in the slice
        """

        if type(curSlice) == slice:
            curIndices = []
            curIndices.append(curSlice.start)
            if curSlice.stop == None:
                # Find the last index
                if dimension == "x" or dimension == "y":
                    dx = collect("dx",\
                                 path=self._path, xguards = self._xguards,\
                                 info=False)
                    dimLen = dx.shape[0]
                if dimension == "y":
                    dy = collect("dy",\
                                 path=self._path, yguards = self._yguards,\
                                 info=False)
                    dimLen = dy.shape[1]
                if dimension == "z":
                    # Subtract 1, as MZ includes the last point
                    dimLen = collect("MZ", path=self._path, info=False) - 1
                if dimension == "t":
                    t = collect("t_array", path=self._path, info=False)
                    dimLen = len(t)
                # Subtract 1 in the end as indices counts from 0
                curIndices.append(dimLen - 1)
            else:
                curIndices.append(curSlice.stop)
        elif curSlice is None:
            curIndices = curSlice
        else:
            curIndices = [curSlice, curSlice]

        # Check for negative indices
        if curIndices is not None:
            for ind in range(len(curIndices)):
                if curIndices[ind] < 0:
                    if dimension == "x" or dimension == "y":
                        dx = collect("dx",\
                                     path=self._path, xguards = self._xguards,\
                                     info=False)
                        dimLen = dx.shape[0]
                    if dimension == "y":
                        dy = collect("dy",\
                                     path=self._path, yguards = self._yguards,\
                                     info=False)
                        dimLen = dy.shape[1]
                    if dimension == "z":
                        # Subtract 1, as MZ includes the last point
                        dimLen = collect("MZ", path=self._path, info=False) - 1
                    if dimension == "t":
                        t   = collect("t_array", path=self._path, info=False)
                        dimLen = len(t)
                    # Subtract 1 in the end as indices counts from 0
                    realInd = dimLen + curIndices[ind] - 1
                    if realInd < 0:
                        message  = ("Index {0} out of range for {1}"
                                    ", as {1} has only {2} elements").\
                            format(curIndices[ind], dimension, dimLen)
                        raise IndexError(message)
                    curIndices[ind] = realInd

        return curIndices
    #}}}
#}}}

#{{{class Plot1D
class Plot1D(Plot):
    """
    Class for plotting the results of the CELMA code in 1D.
    The lines to plot are prepared in the Line class, and the Organizer
    class.

    Inherits from the Plot class.

    Handles:

    * Collection of the variables throug the lines
    * Plotting of the variables
    * Animation of the variables
    """

    #{{{Constructor
    def __init__(self           ,\
                 *args          ,\
                 marker   = None,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Get the proper 1D slice
        3. Sets the marker

        Parameters
        ----------
        *args : positional arguments
            See the constructor of Plot for details.
        marker : str
            The type of marker to be used in the plot.
        **kwargs : keyword arguments
            See the constructor of Plot for details.
        """
        #}}}

        # Call the constructor of the parent class
        super(Plot1D, self).__init__(*args, **kwargs)

        # Check that the indices are properly set
        # Note that this is set after super, as super will check for bad
        # input
        if (type(kwargs["xSlice"]) is slice) and\
           (type(kwargs["ySlice"]) is slice) and\
           (type(kwargs["zSlice"]) is slice):
            message = "3 slices were given, although only 1 is possible"
            raise RuntimeError(message)
        elif (type(kwargs["xSlice"]) == slice and\
              type(kwargs["ySlice"]) == slice) or\
             (type(kwargs["ySlice"]) == slice and\
              type(kwargs["zSlice"]) == slice) or\
             (type(kwargs["zSlice"]) == slice and\
              type(kwargs["xSlice"]) == slice):
            message = "2 slices were given, although only 1 is possible"
            raise ValueError(message)

        # Get the x-axis of the plot
        self._direction = None
        #{{{x-direction
        if type(kwargs["xSlice"]) == slice:
            # Update dict
            self.helper.zTxtDict['value'] =\
                plotNumberFormatter(self.helper.z[kwargs["ySlice"]], None)
            # Set values
            thetaTxt = self.helper.thetaTxtDict["constThetaTxt"].\
                       format(int(self.helper.theta[kwargs["zSlice"]]))
            zTxt     = self.helper.zTxtDict["constZTxt"].\
                       format(self.helper.zTxtDict)
            # Set the label and the title
            self._xAx    = self.helper.rho
            self._xlabel = self.helper.rhoTxtDict["rhoTxtLabel"]
            self._title  = "{}   {}  ".format(thetaTxt, zTxt)

            # Set direction (used in save)
            self._direction = "radial"
        #}}}

        #{{{y-direction
        if type(kwargs["ySlice"]) == slice:
            # Update dict
            self.helper.rhoTxtDict['value'] =\
                plotNumberFormatter(self.helper.rho[kwargs["xSlice"]], None)
            # Set values
            thetaTxt = self.helper.thetaTxtDict["constThetaTxt"].\
                            format(int(self.helper.theta[kwargs["zSlice"]]))
            rhoTxt   = self.helper.rhoTxtDict["constRhoTxt"].\
                            format(self.helper.rhoTxtDict)
            # Set the label and the title
            self._xAx    = self.helper.z
            self._xlabel = self.helper.zTxtDict["zTxtLabel"]
            self._title  = "{}   {}  ".format(rhoTxt, thetaTxt)

            # Set direction (used in save)
            self._direction = "parallel"
        #}}}

        #{{{z-direction
        if type(kwargs["zSlice"]) == slice:

            # Update dicts
            self.helper.rhoTxtDict['value'] =\
                plotNumberFormatter(self.helper.rho[kwargs["xSlice"]], None)
            self.helper.zTxtDict['value'] =\
                plotNumberFormatter(self.helper.z [kwargs["ySlice"]], None)
            # Set values
            rhoTxt = self.helper.rhoTxtDict["constRhoTxt"].\
                            format(self.helper.rhoTxtDict)
            zTxt   = self.helper.zTxtDict["constZTxt"].\
                            format(self.helper.zTxtDict)
            # Set the label and the title
            self._xAx    = r"$\theta$"
            self._xlabel = self.helper.zTxtDict["zTxtLabel"]
            self._title  = "{}   {}   ".format(rhoTxt, zTxt)

            # Set direction (used in save)
            self._direction = "theta"
        #}}}

        if self._direction is None:
            message = ("Improper slicing:\n"
                       "xSlice={}\n"
                       "ySlice={}\n"
                       "zSlice={}\n").format(kwargs["xSlice"],\
                                             kwargs["ySlice"],\
                                             kwargs["zSlice"])
            raise ValueError(message)

        # Set the input data
        self._marker = marker
    #}}}

    #{{{_animFunction
    def _animFunction(self, tInd, orgObj, fig):
        """
        Function which updates the data.

        As blitting is False, there is no reason to return the lines

        Parameters
        ----------
        tInd : int
            The current t index.
        orgObj : Organizer object
            Contains the lines.
        fig : figure
            The figure to plot on.
        """

        # Plot the lines
        for ind, line in enumerate(orgObj.lines):
            # Plot the normal lines
            line.lineObj.set_data(self._xAx, line.field[tInd,:])

            if orgObj.useCombinedPlot:
                # Plot the line in the combined plot
                orgObj.combLineLineObjs[ind].\
                        set_data(self._xAx, line.field[tInd,:])

        # Set the title
        # Update the dictionary
        self.helper.tTxtDict['value'] =\
            plotNumberFormatter(self.helper.t[tInd], None)
        curTimeTxt = self.helper.tTxtDict["tTxt"].format(self.helper.tTxtDict)
        fig.suptitle("{}{}".format(self._title, curTimeTxt))
    #}}}

    #{{{_plotLines
    def _plotLines(self, fig, orgObj, tInd):
        """
        Plots the lines into the combined line plot.

        Parameters
        ----------
        fig : figure
            The figure.
        orgObj : Organizer object
            Contains the lines.
        tInd
            The time index to plot for.
        """

        # Plot the lines, and collect the max and min values
        allMax = []
        allMin = []

        for line in orgObj.lines:
            line.lineObj, =\
                line.ax.plot(self._xAx                     ,\
                             line.field[tInd,:]            ,\
                             marker          = self._marker,\
                             color           = line.color  ,\
                             markeredgecolor = line.color  ,\
                             markerfacecolor = line.color  ,\
                             label           = line.label  ,\
                             )

            # Find the max and the min
            curMax = np.max(line.field)
            curMin = np.min(line.field)

            # Set the y-axis limits
            line.ax.set_ylim(curMin, curMax)

            if orgObj.useCombinedPlot:
                # Store the line object
                orgObj.combLineLineObjs.append(\
                    orgObj.combLine.ax.plot(self._xAx                      ,\
                                             line.field[tInd,:]            ,\
                                             marker          = self._marker,\
                                             color           = line.color  ,\
                                             markeredgecolor = line.color  ,\
                                             markerfacecolor = line.color  ,\
                                             )[0])

                allMax.append(curMax)
                allMin.append(curMin)

            # Decoration
            if line.bottomAx:
                line.ax.set_xlabel(self._xlabel)
            else:
                line.ax.tick_params(labelbottom="off")

            line.ax.legend(loc="upper right", fancybox=True,\
                           framealpha=0.5, numpoints=1)
            # Avoid ticks collision
            line.ax.yaxis.set_major_locator(MaxNLocator(prune="both"))
            line.ax.locator_params(axis="y", tight=True, nbins=6)
            # Avoid silly top value
            line.ax.get_yaxis().get_major_formatter().set_useOffset(False)
            # Use own fuction to deal with ticks
            line.ax.get_yaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))

            line.ax.get_xaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))
            # Set grid
            line.ax.grid(b = True)

        if orgObj.useCombinedPlot:
            # Find the max and the min
            allMax = np.max(allMax)
            allMin = np.min(allMin)

            # Set the y-axis limits
            orgObj.combLine.ax.set_ylim(allMin, allMax)

        # Set the title
        self.helper.tTxtDict['value'] =\
            plotNumberFormatter(self.helper.t[0], None)
        curTimeTxt = self.helper.tTxtDict["tTxt"].format(self.helper.tTxtDict)
        fig.suptitle("{}{}".format(self._title, curTimeTxt))

        # Adjust the subplots
        fig.subplots_adjust(hspace=0, wspace=0.35)
        # Full screen plots
        # http://stackoverflow.com/questions/12439588/how-to-maximize-a-plt-show-window-using-python
        if get_backend() == "QT4Agg":
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
    #}}}

    #{{{collectLine
    def collectLine(self, line):
        """
        Collects the data for one line and reshapes it.

        Parameters
        ----------
        line : Line object
            Line object to set the field to
        """

        if self._subPolAvg:
            # We need to collect the whole field if we would like to do
            # poloidal averages
            try:
                line.field = collect(line.name,\
                                     path    = self._path   ,\
                                     xguards = self._xguards,\
                                     yguards = self._yguards,\
                                     tind    = self._tind   ,\
                                     info    = False)
            except ValueError:
                # Raise an OSError as this is excepted
                raise OSError("Could not collect")

            # If Variable not saved each timestep
            if len(line.field.shape) == 3:
                # Make it a 4d variable
                field      = np.zeros(( len(self.helper.t), *line.field.shape))
                # Copy the field in to each time
                field[:]   = line.field
                line.field = field

            # Subtract the poloidal average, and slice the result
            line.field = (line.field - polAvg(line.field)) \
                    [:,\
                     self._xSlice,\
                     self._ySlice,\
                     self._zSlice,\
                    ]

        else:
            try:
                line.field = collect(line.name,\
                                     path    = self._path   ,\
                                     xguards = self._xguards,\
                                     yguards = self._yguards,\
                                     xind    = self._xind   ,\
                                     yind    = self._yind   ,\
                                     zind    = self._zind   ,\
                                     tind    = self._tind   ,\
                                     info    = False)
            except ValueError:
                # Raise an OSError as this is excepted
                raise OSError("Could not collect")

            # If Variable not saved each timestep
            if len(line.field.shape) == 3:
                # Make it a 4d variable
                field      = np.zeros(( len(self.helper.t), *line.field.shape))
                # Copy the field in to each time
                field[:]   = line.field
                line.field = field

        # Slice in t
        if self._tSlice is not None:
            if self._tSlice.step is not None:
                line.field = line.field[::self._tSlice.step]

        # Flatten the variables except the time dimension
        # -1 => total size divided by product of all other listed dimensions
        line.field = line.field.reshape(line.field.shape[0], -1)
    #}}}

    #{{{plotDriver
    def plotDriver(self, fig, orgObj, savePath = "."):
        """
        Function which drives the plotting.

        Parameters
        ----------
        fig : figure
            The figure to plot on.
        orgObj : Organizer object
            The organization object.
        savePath : str
            Path to save file to.
        """

        # Turn off calculation of physical units if you are not dealing
        # with main fields
        if orgObj.pltName != "mainFields":
            self.convertToPhysical = False

        # Initial plot
        self._plotLines(fig, orgObj, 0)

        if self._savePlot:
            if not os.path.exists(savePath):
                os.makedirs(savePath)
            # Make a saveName by stripping the orgObj's plot name for bad
            # characters
            saveName   = orgObj.pltName.replace("\\", "")
            saveName   = saveName.replace("{", "")
            saveName   = saveName.replace("}", "")
            saveName   = saveName.replace("^", "")
            fileName   = saveName + "-" + self._direction
            fileName = os.path.join(savePath, fileName)

        # Animate if we have more than one frame
        if self._frames > 1:
            anim = animation.FuncAnimation(fig                   ,\
                                           self._animFunction    ,\
                                           fargs  = (orgObj, fig),\
                                           frames = self._frames ,\
                                           blit   = False        ,\
                                           )

            if self._savePlot:
                if self._writer == "ffmpeg":
                    FFMpegWriter = animation.writers['ffmpeg']
                    # * bitrate is set high in order to have ok quality
                    # * fps sets the speed
                    #   http://stackoverflow.com/questions/22010586/matplotlib-animation-duration
                    # * codec is by default mpeg4, but as this creates large
                    #   files. h264 is preferred.
                    # * For installation, see
                    #   https://github.com/loeiten/usingLinux/blob/master/installationProcedures/ffmpeg.md
                    self._writer = FFMpegWriter(bitrate = self._bitrate,\
                                                fps     = self._fps    ,\
                                                codec   = self._codec)
                # Save the animation
                anim.save(fileName + self._animExtension   ,\
                          writer = self._writer            ,\
                          savefig_kwargs = {"pad_inches":0},\
                          )
                print("Saved to {}{}".format(fileName, self._animExtension))
        else:
            if self._savePlot:
                # Save the figure
                fig.savefig("{}.{}".format(fileName, self._extension),\
                            transparent = False  ,\
                            bbox_inches = "tight",\
                            pad_inches  = 0      ,\
                            )
                print("Saved to {}.{}".format(fileName, self._extension))

        if self._showPlot:
            fig.show()

        plt.close(fig)
    #}}}
#}}}

#{{{class Plot2D
class Plot2D(Plot):
    """
    Class for plotting the results of the CELMA code in 2D.
    Inherits from the Plot class

    Handles:

    * Collection of the variables
    * Plotting of the variables
    * Animation of the variables
    """

    #{{{Constructor
    def __init__(self                            ,\
                 path                            ,\
                 varName                         ,\
                 var               = None        ,\
                 maxGradRhoFolder  = None        ,\
                 xguards           = False       ,\
                 yguards           = False       ,\
                 varMax            = None        ,\
                 varMin            = None        ,\
                 varyMaxMin        = False       ,\
                 axisEqualParallel = True        ,\
                 mode              = "perpAndPar",\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Builds the mesh from the CylinderMesh class
        2. Either:
            1. Collects the variable
            2. Calculates the variable form varFunc if set
        3. Calculates the physical units if convertToPhysical is set
        4. Sets the lines which shows how the data is sliced
        5. Initializes the plot

        Parameters
        ----------
        path : str
            Path to collect from
        varName : str
            Name of the field which is going to be collected (if var is
            not given).
        var : [None | array]
            The variable to plot.
        maxGradRhoFolder : [None | str]
            If this is set, the xSlice will be replaced by the index of
            the largest gradient in rho direction.
        varMax : [None | float]
            Setting a hard upper limit z-axis in the plot.
        varMin : [None | float]
            Setting a hard lower limit z-axis in the plot.
        varyMaxMin : bool
            Whether or not the limits of the z-axis should be
            set to the max/min of the current timestep or not.
        axisEqualParallel : bool
            Whether or not the parallel plot should be plotted with axis
            equal or not.
        mode : ["perpAndPar" | "perpAndPol" | "perp" | "par" | "pol"]
            The mode to plot.
        **kwargs : keyword arguments
            See the constructor of Plot for details.
        """
        #}}}

        # Check that mode is set correctly:
        implementedModes = ("perpAndPar",  "perpAndPol", "perp", "par", "pol")
        found = False
        for checkMode in implementedModes:
            if mode == checkMode:
                self._mode = mode.lower()
                found = True
                break

        if not found:
            message = "mode '{}' not implemented".format(mode)
            raise NotImplementedError(message)

        # Call the constructor of the parent class
        super(Plot2D, self).__init__(path   ,\
                                     xguards,\
                                     yguards,\
                                     **kwargs)

        # Check that the indices are properly set
        if (self._xSlice == slice(0,None)) and\
           (self._ySlice == slice(0,None)) and\
           (self._zSlice == slice(0,None)):
            message = "3 slices were given, although only 2 is possible"
            raise ValueError(message)

        #{{{Reset xSlice through maxGradRhoFolder
        if maxGradRhoFolder:
            # Collect the profile of the variable
            if varName == "n":
                curVarName = "lnN"
            else:
                curVarName = varName
            # Need the x-guards as derivatives will be taken (subtracted
            # in findLargestRadialGrad
            tmpVar = collect(varname = curVarName                  ,\
                             path    = maxGradRhoFolder            ,\
                             xguards = True                        ,\
                             yguards = False                       ,\
                             info    = False                       ,\
                             yind    = (self._ySlice, self._ySlice),\
                             zind    = (self._zSlice, self._zSlice),\
                             )
            if varName == "n":
                tmpVar = np.exp(tmpVar)
            # Collect variables needed for findLargestRadialGrad
            dx  = collect(varname = "dx"                        ,\
                          path    = path                        ,\
                          xguards = True                        ,\
                          yguards = False                       ,\
                          info    = False                       ,\
                          yind    = (self._ySlice, self._ySlice),\
                          )
            MXG = collect(varname = "MXG", path = path, info=False)
            # Find the max gradient of the variable (subtracts the guard cells)
            _, self._xSlice = findLargestRadialGrad(tmpVar, dx, MXG)
            # Update _xind
            self._xind = self._getIndices(self._xSlice, "x")
        #}}}

        # Make it possible to filter warnings (Ex: no variation in the data)
        warnings.filterwarnings("error")

        # Set member data from the index
        self._varyMaxMin        = varyMaxMin
        self._varMax            = varMax
        self._varMin            = varMin
        self._axisEqualParallel = axisEqualParallel

        # Set additional plot properties
        self._latexSize = 35
        self._nCont     = 100
        self._pltName   = None

        # Create a CylinderMesh object
        self._cyl = CylinderMesh(self.helper.rho,\
                                 self.helper.theta,\
                                 self.helper.z,\
                                 xguards)

        # Collect the full variable
        # Stored as an ndarray with the indices [t,x,y,z] (=[t,rho,z,theta])
        if var is None:
            self._variable = collect(varName             ,\
                                     path    = path      ,\
                                     yguards = yguards   ,\
                                     xguards = xguards   ,\
                                     tind    = self._tind,\
                                     info    = False     ,\
                                     )
        else:
            self._variable = var

        # Slice in t
        if self._tSlice is not None:
            if self._tSlice.step is not None:
                self._variable = self._variable[::self._tSlice.step]

        if self._subPolAvg:
            self._variable = self._variable - polAvg(self._variable)

        # Add the last theta slice
        self._variable =\
                self._cyl.addLastThetaSlice(self._variable, len(self.helper.t))

        self._variable, self._normalization, self._units =\
                self.helper.physicalUnitsConverter(self._variable, varName)

        if xguards:
            # Remove the inner ghost points from the variable
            self._variable = np.delete(self._variable, (0), axis=1)

        # Get the max and the min so that we can keep the color coding correct
        if self._varMax == None:
            self._varMax = np.max(self._variable)
        if self._varMin == None:
            self._varMin = np.min(self._variable)
        # Diverging colormap for fluctuations
        if self._subPolAvg:
            self._varMax = np.max([np.abs(self._varMax), np.abs(self._varMin)])
            self._varMin = - self._varMax

        # We need to manually sepcify the levels in order to have a
        # fixed color bar
        self._levels = np.linspace(self._varMin  ,\
                                   self._varMax  ,\
                                   self._nCont   ,\
                                   endpoint = True)

        # Then theta index corresponding to pi
        piInd = round(self._variable.shape[3]/2)

        # Calculate the theta in degrees
        self._thetaRad = self.helper.dz*self._zind[-1]
        self._thetaDeg = self._thetaRad*(180/np.pi)

        # Get the Z values of the X, Y, Z plots
        # We subscript the last index of self._@ind, as this is given as
        # a range in the Plot constructor
        if "pol" in self._mode:
            self._Z_ZT = self._variable[:, self._xind[-1],: , :]
        if "perp" in self._mode:
            self._Z_RT = self._variable[:, :, self._yind[-1], :]
        if "par" in self._mode:
            self._Z_RZ = self._variable[:, :, :, self._zind[-1]]
            # Get the Z value in the RZ plane which is pi from the current index
            if self._zind[-1] > piInd:
                self._Z_RZ_P_PI = self._variable[:, :, :, self._zind[-1] - piInd]
            else:
                self._Z_RZ_P_PI = self._variable[:, :, :, self._zind[-1] + piInd]

        if self._mode == "perpAndPar".lower()\
           or self._mode == "perpAndPol".lower():
            self._createLines()

        # Create the figure and axis
        if self._mode == "perpAndPar".lower() or self._mode == "perpAndPol".lower():
            pltSize = (30,15)
        else:
            pltSize = (20,15)
        self._fig    = plt.figure(figsize = pltSize)
        #{{{if self._mode == "perpAndPar".lower()
        if self._mode == "perpAndPar".lower():
            gs           = GridSpec(1, 3, width_ratios=(20, 20, 1))
            self._perpAx = self._fig.add_subplot(gs[0])
            self._parAx  = self._fig.add_subplot(gs[1])
            self._cBarAx = self._fig.add_subplot(gs[2])
            self._fig.subplots_adjust(wspace=0.25)
            self._parAx.grid(True)
            self._perpAx.grid(True)
        #}}}
        #{{{elif self._mode == "perpAndPol".lower()
        elif self._mode == "perpAndPol".lower():
            gs           = GridSpec(1, 3, width_ratios=(20, 20, 1))
            self._perpAx = self._fig.add_subplot(gs[0])
            self._polAx  = self._fig.add_subplot(gs[1])
            self._cBarAx = self._fig.add_subplot(gs[2])
            self._fig.subplots_adjust(wspace=0.25)
            self._perpAx.grid(True)
            self._polAx.grid(True)
        #}}}
        #{{{elif self._mode == "perp"
        elif self._mode == "perp":
            self._perpAx = self._fig.add_subplot(111)
            self._cBarAx = make_axes_locatable(self._perpAx).\
                           append_axes('right', '5%', '5%')
            self._perpAx.grid(True)
        #}}}
        #{{{elif self._mode == "par"
        elif self._mode == "par":
            self._parAx = self._fig.add_subplot(111)
            self._cBarAx = make_axes_locatable(self._parAx).\
                           append_axes('right', '5%', '5%')
            self._parAx.grid(True)
        #}}}
        #{{{elif self._mode == "pol"
        elif self._mode == "pol":
            self._polAx = self._fig.add_subplot(111)
            self._cBarAx = make_axes_locatable(self._polAx).\
                           append_axes('right', '5%', '5%')
            self._polAx.grid(True)
        #}}}

        # Set memberdatas altered in _plot2D
        self._cbarPlane = None
        #}}}

    #{{{_createLines
    def _createLines(self):
        """ Set the lines which shows where the data is sliced"""

        if self._mode == "perpAndPar".lower():
            # The slice lines we are plotting
            rhoStart = self.helper.rho[0]
            rhoEnd   = self.helper.rho[-1]

            # Calculate the numerical value of the theta angle and the z value
            thetaRad = self._thetaRad
            thetaPPi = thetaRad + np.pi
            zVal     = self.helper.z[self._ySlice]

            # Set coordinates for the lines which indicate how the data is
            # sliced
            # Organized in start and stop pairs
            # We need two lines due to the center of the cylinder
            self._RTLine1XVals =\
                    (rhoStart*np.cos(thetaRad), rhoEnd*np.cos(thetaRad))
            self._RTLine1YVals =\
                    (rhoStart*np.sin(thetaRad), rhoEnd*np.sin(thetaRad))
            self._RTLine2XVals =\
                    (rhoStart*np.cos(thetaPPi), rhoEnd*np.cos(thetaPPi))
            self._RTLine2YVals =\
                    (rhoStart*np.sin(thetaPPi), rhoEnd*np.sin(thetaPPi))
            self._RZLine1XVals =\
                    (-rhoEnd                  , -rhoStart              )
            self._RZLine1YVals =\
                    (zVal                     , zVal                   )
            self._RZLine2XVals =\
                    (rhoStart                 , rhoEnd                 )
            self._RZLine2YVals =\
                    (zVal                     , zVal                   )
        elif self._mode == "perpAndPol".lower():
            # Get the radius of the circle
            rhoVal = self.helper.rho[self._xSlice]

            # Create the circle
            self._circle = plt.Circle((0, 0)         ,\
                                      rhoVal         ,\
                                      fill = False   ,\
                                      linewidth = 1.0,\
                                      linestyle= '--',\
                                      color='k')

            # Get the z value
            zVal = self.helper.z[self._ySlice]

            # Set coordinates for the lines which indicate how the data is
            # sliced
            # Organized in start and stop pairs
            # We only need one line to indicate the z value on the
            # poloidal plot
            self._ZTLineXVals=(0   , 2.0*np.pi)
            self._ZTLineYVals=(zVal, zVal     )
    #}}}

    #{{{_plot2D
    def _plot2D(self, tInd):
        #{{{docstring
        """
        Performs the actual plotting

        Parameters
        ----------
        tInd : int
            The index to plot for.
        """
        #}}}
        # Check that levels are rising
        if not(self._levels is None):
            if len(self._levels) > 1 and np.amin(np.diff(self._levels)) <= 0.0:
                self._levels = None

        if self._varyMaxMin:
            # Allocate the max and min lists
            maxList = []
            minList = []
        #{{{if "pol" in self._mode
        if "pol" in self._mode:
            Z_ZT = self._Z_ZT[tInd, :, :]
            if self._varyMaxMin:
                maxList.append(np.max(Z_ZT))
                minList.append(np.min(Z_ZT))
        #}}}
        #{{{if "perp" in self._mode
        if "perp" in self._mode:
            Z_RT = self._Z_RT[tInd, :, :]
            if self._varyMaxMin:
                maxList.append(np.max(Z_RT))
                minList.append(np.min(Z_RT))
        #}}}
        #{{{if "par" in self._mode
        if "par" in self._mode:
            Z_RZ      = self._Z_RZ     [tInd, :, :]
            Z_RZ_P_PI = self._Z_RZ_P_PI[tInd, :, :]
            if self._varyMaxMin:
                maxList.append(np.max(Z_RZ))
                maxList.append(np.max(Z_RZ_P_PI))
                minList.append(np.min(Z_RZ))
                minList.append(np.min(Z_RZ_P_PI))
        #}}}

        # If we want the max and min to vary
        if self._varyMaxMin and tInd:
            # Update the max and min
            self._varMax = np.max(maxList)
            self._varMin = np.min(minList)

            # Diverging colormap for fluctuations
            if self._subPolAvg:
                self._varMax =\
                        np.max([np.abs(self._varMax), np.abs(self._varMin)])
                self._varMin = - self._varMax

            # Update the levels
            levels = np.linspace(self._varMin   ,\
                                 self._varMax   ,\
                                 self._nCont    ,\
                                 endpoint = True,\
                                 )
            # Update the levels just if there is any difference
            if np.amin(np.diff(levels)) > 0.0:
                self._levels = levels

        # Specify repeated kwargs of contourf
        # NOTE: It doens't make sense to use functools.partial here as
        #       contourf is a memberfunction of ax
        # NOTE: zorder sets the rasterization
        #       http://stackoverflow.com/questions/37020842/reducing-size-of-vectorized-contourplot
        cfKwargs = {"cmap"   : self._cmap  ,\
                    "vmax"   : self._varMax,\
                    "vmin"   : self._varMin,\
                    "levels" : self._levels,\
                    "zorder" : -20         ,\
                   }

        #{{{ Plot, set labels and draw grids
        #{{{if "pol" in self._mode
        if "pol" in self._mode:
            # Clear previous axis
            self._polAx.cla()
            # Plot the poloidal plane
            polPlane = self._polAx.\
                contourf(\
                self._cyl.X_ZT, self._cyl.Y_ZT, Z_ZT.transpose(), **cfKwargs)
            # Set rasterization order
            self._polAx.set_rasterization_zorder(-10)
            # Draw the grids
            self._polAx.grid(b=True)
            # Set x and y labels
            self._polAx.set_xlabel(r"$\theta$", fontsize = self._latexSize)
            self._polAx.set_ylabel(self.helper.zTxtDict["zTxtLabel"],\
                                   fontsize = self._latexSize)
        #}}}
        #{{{if "perp" in self._mode
        if "perp" in self._mode:
            # Clear previous axis
            self._perpAx.cla()
            # Plot the perpendicular plane
            perpPlane = self._perpAx.\
                    contourf(self._cyl.X_RT, self._cyl.Y_RT, Z_RT, **cfKwargs)
            # Set rasterization order
            self._perpAx.set_rasterization_zorder(-10)
            # Draw the grids
            self._perpAx.grid(b=True)
            # Set x and y labels
            self._perpAx.set_xlabel(self.helper.rhoTxtDict["rhoTxtLabel"],\
                                 fontsize = self._latexSize)
            self._perpAx.set_ylabel(self.helper.rhoTxtDict["rhoTxtLabel"],\
                                 fontsize = self._latexSize)
        #}}}
        #{{{if "par" in self._mode
        if "par" in self._mode:
            # Clear previous axis
            self._parAx.cla()
            # Plot the parallel plane
            parPlane  = self._parAx.\
            contourf(self._cyl.X_RZ, self._cyl.Y_RZ, Z_RZ, **cfKwargs)
            # parPlaneNeg, not used
            self._parAx.\
            contourf(self._cyl.X_RZ_NEG, self._cyl.Y_RZ, Z_RZ_P_PI, **cfKwargs)
            # Set rasterization order
            self._parAx.set_rasterization_zorder(-10)
            # Draw the grids
            self._parAx.grid(b=True)
            # Set x and y labels
            self._parAx.set_xlabel(self.helper.rhoTxtDict["rhoTxtLabel"],\
                                 fontsize = self._latexSize)
            self._parAx.set_ylabel(self.helper.zTxtDict["zTxtLabel"],\
                                fontsize = self._latexSize)
        #}}}
        #}}}

        #{{{ Draw lines
        #{{{if self._mode == "perpAndPar"
        if self._mode == "perpAndPar".lower():
            # Lines needs only to be plotted once
            # Par line 1
            self._parAx.plot(self._RZLine1XVals,\
                           self._RZLine1YVals,\
                           "--k"             ,\
                           linewidth = 1     ,\
                           )
            # Par line 2
            self._parAx.plot(self._RZLine2XVals,\
                           self._RZLine2YVals,\
                           "--k"             ,\
                           linewidth = 1     ,\
                           )
            # Perp line 1
            self._perpAx.plot(self._RTLine1XVals,\
                           self._RTLine1YVals,\
                           "--k"             ,\
                           linewidth = 1     ,\
                           )
            # Perp line 2
            self._perpAx.plot(self._RTLine2XVals,\
                           self._RTLine2YVals,\
                           "--k"             ,\
                           linewidth = 1     ,\
                           )

        #}}}
        #{{{elif self._mode == "perpAndPol"
        elif self._mode == "perpAndPol".lower():
            # Lines needs only to be plotted once
            # Circle
            self._perpAx.add_artist(self._circle)

            # Pol line
            self._polAx.plot(self._ZTLineXVals,\
                             self._ZTLineYVals,\
                             "--k"            ,\
                             linewidth = 1    ,\
                             )
        #}}}
        #}}}


        # Title preparation
        self.helper.rhoTxtDict["value"] =\
                plotNumberFormatter(self.helper.rho[self._xSlice], None)
        self.helper.zTxtDict["value"] =\
                plotNumberFormatter(self.helper.z[self._ySlice], None)
        self.helper.tTxtDict["value"] =\
                plotNumberFormatter(self.helper.t[tInd], None, precision=4)

        # Titles
        polTitle =\
            self.helper.rhoTxtDict["constRhoTxt"].format(self.helper.rhoTxtDict)
        perpTitle =\
            self.helper.zTxtDict["constZTxt"].format(self.helper.zTxtDict)
        parTitle =\
        self.helper.thetaTxtDict["constThetaTxt"].format(int(self._thetaDeg))
        timeTitle = self.helper.tTxtDict["tTxt"].format(self.helper.tTxtDict)

        # Specify repeated kwargs of txt
        txtKwargs = { "horizontalalignment" : "center"       ,\
                      "verticalalignment"   : "center"       ,\
                      "fontsize"            : self._latexSize,\
                    }
        #{{{ Set the titles
        #{{{if self._mode == "perpAndPar".lower()
        if self._mode == "perpAndPar".lower():
            # Perp text
            self._perpAx.text(0.5, 1.05, perpTitle,\
                         transform = self._perpAx.transAxes,\
                         **txtKwargs)

            # Par text
            self._parAx.text(0.5, 1.05, parTitle,\
                         transform = self._parAx.transAxes,\
                         **txtKwargs)

            # Title mid
            # Text for the figure. Could append this to the figure itself,
            # but it seems to be easier to just add it to an axis due to
            # animation
            self._figTxt = self._perpAx.text(1.10, 1.05, timeTitle,\
                                          transform = self._perpAx.transAxes,\
                                          **txtKwargs)
        #}}}
        #{{{elif self._mode == "perpAndPol".lower()
        elif self._mode == "perpAndPol".lower():
            # Perp text
            self._perpAx.text(0.5, 1.05, perpTitle,\
                         transform = self._perpAx.transAxes,\
                         **txtKwargs)

            # Pol text
            self._polAx.text(0.5, 1.05, polTitle,\
                         transform = self._polAx.transAxes,\
                         **txtKwargs)

            # Title mid
            # Text for the figure. Could append this to the figure itself,
            # but it seems to be easier to just add it to an axis due to
            # animation
            self._figTxt = self._perpAx.text(1.10, 1.05, timeTitle,\
                                          transform = self._perpAx.transAxes,\
                                          **txtKwargs)
        #}}}
        #{{{elif self._mode == "pol"
        elif self._mode == "pol":
            self._polTxt = self._polAx.text(0.5, 1.05,\
                                         "{}$,$ {}".format(polTitle, timeTitle),\
                                         transform = self._polAx.transAxes,\
                                         **txtKwargs)
        #}}}
        #{{{elif self._mode == "perp"
        elif self._mode == "perp":
            self._perpTxt = self._perpAx.text(0.5, 1.05,\
                                         "{}$,$ {}".format(perpTitle, timeTitle),\
                                         transform = self._perpAx.transAxes,\
                                         **txtKwargs)
            self._txtSet = True
        #}}}
        #{{{elif self._mode == "par"
        elif self._mode == "par":
            self._parTxt = self._parAx.text(0.5, 1.05,\
                                         "{}$,$ {}".format(parTitle, timeTitle),\
                                         transform = self._parAx.transAxes,\
                                         **txtKwargs)
            self._txtSet = True
        #}}}
        #}}}

        #{{{ Format axes and set equal
        #{{{if "pol" in self._mode:
        if "pol" in self._mode:
            self._polAx.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
            self._polAx.set_xticklabels(\
                    [r"$0$", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
            self._polAx.get_yaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))
        #}}}
        #{{{if "perp" in self._mode:
        if "perp" in self._mode:
            self._perpAx.get_xaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))
            self._perpAx.get_yaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))
            self._perpAx.axis("equal")
        #}}}
        #{{{if "par" in self._mode:
        if "par" in self._mode:
            self._parAx.get_xaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))
            self._parAx.get_yaxis().set_major_formatter(\
                FuncFormatter(plotNumberFormatter))
            if self._axisEqualParallel:
                self._parAx.axis("equal")
        #}}}
        #}}}

        #{{{ Set self._cbarPlane
        if "perp" in self._mode:
            self._cbarPlane = perpPlane
        elif "par" in self._mode:
            self._cbarPlane = parPlane
        elif "pol" in self._mode:
            self._cbarPlane = polPlane
        #}}}

        # Set the colorbar
        # Clear the axis
        # http://stackoverflow.com/questions/39472017/how-to-animate-the-colorbar-in-matplotlib/39596853
        self._cBarAx.cla()
        self._fig.colorbar(self._cbarPlane, cax = self._cBarAx,\
                           format = FuncFormatter(plotNumberFormatter))
    #}}}

    #{{{plotDriver
    def plotDriver(self, pltName, savePath = "."):
        """
        Function which drived the plotting.

        Parameters
        ----------
        pltName : str
            Name of the plot written in LaTeX format, but without the $.
        savePath : str
            Path to save file to.
        """

        self._pltName = pltName

        # Initial plot
        self._plot2D(0)

        if self._mode == "perpAndPar".lower() or\
           self._mode == "perpAndPol".lower():
            # Need to specify rect in order to have top text
            self._fig.tight_layout(w_pad = 2.5, rect=(0,0,1,0.97))

        if self._savePlot:
            # Make dir if not exists
            if not os.path.exists(savePath):
                os.makedirs(savePath)
            # Make a saveName by stripping the orgObj's plot name for bad
            # characters
            saveName = pltName.replace("\\", "")
            saveName = saveName.replace("{", "")
            saveName = saveName.replace("}", "")
            saveName = saveName.replace("^", "")
            fileName = "{}-{}-{}".format(saveName, self._mode, "2D")
            fileName = os.path.join(savePath, fileName)

        # Animate if we have more than one frame
        if self._frames > 1:

            # Animate
            anim = animation.FuncAnimation(self._fig            ,\
                                           self._plot2D         ,\
                                           frames = self._frames,\
                                           blit   = False       ,\
                                           )

            if self._savePlot:
                if self._writer == "ffmpeg":
                    FFMpegWriter = animation.writers['ffmpeg']
                    # * bitrate is set high in order to have ok quality
                    # * fps sets the speed
                    #   http://stackoverflow.com/questions/22010586/matplotlib-animation-duration
                    # * codec is by default mpeg4, but as this creates large
                    #   files. h264 is preferred.
                    # * For installation, see
                    #   https://github.com/loeiten/usingLinux/blob/master/installationProcedures/ffmpeg.md
                    self._writer = FFMpegWriter(bitrate = self._bitrate,\
                                                fps     = self._fps    ,\
                                                codec   = self._codec)
                # Save the animation
                anim.save(fileName + self._animExtension,\
                          writer = self._writer         ,\
                          )
                print("Saved to {}{}".format(fileName, self._animExtension))
        else:
            if self._savePlot:
                # Save the figure
                self._fig.savefig("{}.{}".format(fileName, self._extension),\
                                  transparent = False  ,\
                                  bbox_inches = "tight",\
                                  pad_inches  = 0      ,\
                                  )
                print("Saved to {}.{}".format(fileName, self._extension))

        if self._showPlot:
            self._fig.show()

        plt.close(self._fig)
    #}}}
#}}}
