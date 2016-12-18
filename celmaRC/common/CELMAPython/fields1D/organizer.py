#!/usr/bin/env python

from .line import Line
from ..plotHelpers import seqCMap2
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from boututils.datafile import DataFile
import os


"""
Contains the organizer class
"""




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
from ..plotHelpers import safeCollect
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


# FIXME: This part is in plot2DSuperClass
        # Set the bitrates, fps and codec for ffmpeg (currently magic numbers)
        self._bitrate = -1
        self._fps     = 10
        self._codec   = "h264"

        # Set member data from input
        self._xguards    = xguards
        self._yguards    = yguards

# FIXME: This part is in plot2DSuperClass
        self._savePlot   = savePlot
        self._saveFolder = saveFolder
        self._path       = path
        self._fluctuation  = fluctuation
        self._extension  = extension
        # Set colormap
        if self._fluctuation:
            self._cmap = divCMap
        else:
            self._cmap = seqCMap

# FIXME: These are in collectAndCalc2D
        # Get proper indices
        self._xind = self.slicesToIndices(xSlice, "x")
        self._yind = self.slicesToIndices(ySlice, "y")
        self._zind = self.slicesToIndices(zSlice, "z")
        self._tind = self.slicesToIndices(tSlice, "t")


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

        # Get the x-axis of the plot
        self._direction = None
# FIXME: Here obtain the plot label and text
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

        if self._fluctuation:
            # We need to collect the whole field if we would like to do
            # poloidal averages
            try:
                line.field = safeCollect(line.name,\
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
                line.field = safeCollect(line.name,\
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

orgObj??
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




# FIXME: Is org plot or collect?
#{{{class Organizer
class Organizer(object):
    """
    Class which organizes several lines in a 1D plot.

    This class is responsible for

    * Organization of the lines
    * Make all axes
    * Setting proper names on figures and files
    * Making a combined line if useCombinedPlot is true
    """

    #{{{Constructor
    def __init__(self                    ,\
                 pltName                 ,\
                 combLineName    = None  ,\
                 cols            = 2     ,\
                 useCombinedPlot = False ,\
                 forceCombined   = True  ,\
                 path            = "data",\
                 ):
        """
        The constructor initializes the sequence of lines

        Parameters
        ----------
        pltName : str
            Name of the plot written in LaTeX format, but without the $
        combLineName : [None | str]
            If the string is set, then the organizer will try to collect
            ddt(combLineName) in makeCombinedLine.
        cols : int
            The total number of columns to be used in the plot
        useCombinedPlot : bool
            Toggles if a plot of combined lines are to be plotted
        forceCombined : bool
            Will still plot combined lines, even if lines are not found
        path : str
            Path to collect from
        """

        # Set member data from input
        self._cols           = cols
        self.useCombinedPlot = useCombinedPlot
        self.pltName         = pltName
        self.combLineName    = combLineName
        self._forceCombined   = forceCombined

        # Initialize non-input members
        self._pltSize         = (18,12)
        self.combLine         = None
        self.combLineLineObjs = []
        self.lines            = []
        self.extraLines       = {}
        self.axes             = []


    #}}}

    #{{{pltPrepare
    def pltPrepare(self):
        """
        Prepares the lines in a plot for plotting.
        Call this function before calling collect.

        1. Check that a line is collectable.
        2. Add extra lines (combination of fields etc.) if any.
        3. Check if any plot pos have been given.
           If yes, lines will be rearranged.
        4. Set the color of each plot.
        5. Finds the bottom axes.
        6. Creates the figure and axes.
        7. Returns the figure.
        """






        # Check that lines are collectable
#{{{

        # Variables collectable in the dump file
        dataFileVars = DataFile(os.path.join(path,"BOUT.dmp.0.nc")).list()
        # Make everything lowercase in order to easen comparison
        self._dataFileVars = tuple(el.lower() for el in dataFileVars)

        notFound = []
        for line in self.lines:
            lowerCaseLin = line.name.lower()
            if lowerCaseLin not in self._dataFileVars:
                if self._forceCombined:
                    message = ("{0}!!! Warning: {1} could not be found. "
                               "A combined line will still be plotted as "
                               "forceCombined is set to True.{0}")
                else:
                    message = "{0}!!! Warning: {1} could not be found{0}"
                print(message.format("\n"*2, line.name))
                notFound.append(line)
        for missing in notFound:
            self.lines.remove(missing)
            if not(self._forceCombined):
                self.useCombinedPlot = False
#}}}







        # Append lines with eventual extra lines
        # The extra lines are treated specially since they are not
        # collectable
        nExtraLines = len(self.extraLines.keys())
        if nExtraLines > 0:
            for key in self.extraLines.keys():
                # Insert the extraLines into self.lines
                self.lines.append(self.extraLines[key])

        # Organize the lines
        newLines = [None]*len(self.lines)
        # Make a copy as we are going to pop the list
        self.allFieldsPresent = True
        for line in self.lines:
            if line.plotPos:
                # Get index in self line
                index = self.lines.index(line)
                try:
                    newLines[line.plotPos] = line
                except IndexError:
                    # If not all the fields are saved
                    print(("WARNING: Could not properly position\n\n"))
                    self.allFieldsPresent = False
                    break
        if len(self.lines) > 0 and self.allFieldsPresent:
            # Get the free indices in newLines
            indices = tuple(i for i, el in enumerate(newLines) if el is None)
            for index, line in zip(indices, self.lines):
                newLines[index] = line

        if self.allFieldsPresent:
            # Reassign
            self.lines = newLines
        else:
            print(("Removing extra lines as positioning failed"))
            for key in self.extraLines.keys():
                # Insert the extraLines into self.lines
                self.lines.remove(self.extraLines[key])

        # Set the colors
        colorSpace = np.arange(len(self.lines))
        colors = seqCMap2(np.linspace(0, 1, len(colorSpace)))

        for lineNr, line in enumerate(self.lines):
            line.color = colors[lineNr]

        # If a combined line is to be plotted
        if self.useCombinedPlot:
            # Make a line object
            self.combLine = Line(name  = "ddt({})".format(self.combLineName),\
                                 label = r"\partial_t " + self.pltName ,\
                                 )
            # Make the lastline black, and append it to the lines
            self.combLine.color = "k"
            self.lines.append(self.combLine)

        # Calculate the number of rows
        rows = int(np.ceil(len(self.lines)/self._cols))

        # Create the figure
        fig = plt.figure(figsize = self._pltSize)
        gs  = GridSpec(rows, self._cols)

        # Make the axes
        for lineNr, line in enumerate(self.lines):
            if lineNr == 0:
                # Need an initial line
                line.ax = fig.add_subplot(gs[lineNr])
                firstAx = line.ax
            else:
                line.ax = fig.add_subplot(gs[lineNr], sharex=firstAx)

        for col in range(1, self._cols+1):
            try:
                self.lines[-col].bottomAx = True
            except IndexError:
                # Only one column
                if self.useCombinedPlot:
                    print("WARNING: Only combLine found. Will not plot!!!\n\n")
                    plt.close(fig)
                    return None
                break

        # Pop the combined line in order not to collect it
        if self.useCombinedPlot:
            self.combLine = self.lines.pop()

        # Pop the extra lines in order not to collect them
        if nExtraLines > 0 and self.allFieldsPresent:
            for key in self.extraLines.keys():
                # Update the extraLines, and remove them from
                # self.lines, in order not to collect them
                ind = self.lines.index(self.extraLines[key])
                self.extraLines[key] = self.lines.pop(ind)

        fig.canvas.set_window_title(self.pltName)

        return fig
    #}}}

    #{{{makeCombinedLine
    def makeCombinedLine(self, plotter):
        """
        Makes a combined line.
        To be called after all other lines are collected.

        The routine will first try to collect ddt of the labelName. If
        this is not available, a sum of the lines is used instead.

        Parameters
        ----------
        plotter : plotter object
            Containing the collectLine function
        """

        try:
            plotter.collectLine(self.combLine)
            print("ddt was collected")
        except OSError:
            # OSError is thrown if filed is not found

            # Initialize the field
            self.combLine.field = np.zeros(self.lines[0].field.shape)

            for line in self.lines:
                self.combLine.field += line.field

            print("ddt was obtained by summation")

        # Re-add the combLine to the list
        self.lines.append(self.combLine)
    #}}}
#}}}















#!/usr/bin/env python

"""
Contains the Line class.
"""

#{{{class Line
class Line(object):
    """
    Class containing the data for the line which are going to be
    plotted.

    Provides the data which is connected to a line in a subplot.
    Specifically it contains information about:

    * name     - The name used for collecting the data from a simulation
    * label    - The label which is going to be used for in the legend
    * plotPos  - Index of the plot number (if any)
    * ax       - The axis object for the line
    * bottomAx - If the axis is the lowest
    * lineObj  - The lineobject of the plot
    * color    - The color of the line
    * field    - The data used to plot the line
    * dim      - The dimension of the data
    """

    #{{{Constructor
    def __init__(self          ,\
                 name          ,\
                 label         ,\
                 plotPos = None,\
                ):
        #{{{docstring
        """
        The constructor sets the member data

        Parameters
        ----------
        name : str
            Name which is going to be collected
        label : str
            Label which is going to be used in the plot
            NOTE: "$" will be added around this string
            NOTE: Remember to use raw strings
        plotPos : int
            If there is a preferred position of the plot.
            Given as an index number.
            NOTE: The user have to take care so that two plots
            does not share the same number.
        """
        #}}}

        # Set member data from input
        self.name    = name
        self.label   = r"$" + label + r"$"
        self.plotPos = plotPos

        # Initialize non-input members
        self.ax       = None
        self.bottomAx = False
        self.lineObj  = None
        self.color    = None
        self.field    = None
        self.dim      = None
    #}}}
#}}}
