#!/usr/bin/env python

"""
Contains classes for plotting the fields
"""

from .. import titleSize
from ..statsAndSignals import polAvg
from .cylinderMesh import CylinderMesh
import scipy.constants as cst
from matplotlib import get_backend
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from boutdata import collect
from boututils.options import BOUTOptions
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
    def __init__(self                      ,\
                 path                      ,\
                 xguards    = False        ,\
                 yguards    = False        ,\
                 xSlice     = slice(0,None),\
                 ySlice     = slice(0,None),\
                 zSlice     = slice(0,None),\
                 tSlice     = None         ,\
                 physicalU  = False        ,\
                 subPolAvg  = False        ,\
                 showPlot   = False        ,\
                 savePlot   = True         ,\
                 saveFolder = None         ,\
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
        physicalU : bool
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
        """
        #}}}

        # Set member data from input
        self._path       = path
        self._xguards    = xguards
        self._yguards    = yguards
        self._showPlot   = showPlot
        self._savePlot   = savePlot
        self._saveFolder = saveFolder
        # Public as used in the driver
        self.physicalU   = physicalU

        # Get the coordinates
        #{{{rho
        dx = collect("dx"             ,\
                     path    = path   ,\
                     xguards = xguards,\
                     yguards = yguards,\
                     info    = False)
        MXG = collect("MXG"            ,\
                      path    = path   ,\
                      xguards = xguards,\
                      yguards = yguards,\
                      info    = False)

        nPoints = dx.shape[0]
        dx      = dx[0,0]

        if xguards:
            innerPoints = nPoints - 2*MXG
        else:
            innerPoints = nPoints

        # By default there is no offset in the cylinder
        # For comparision with other codes, an offset option is set
        # Read the input file
        myOpts = BOUTOptions(path)
        # Read in geom offset
        try:
            offset = eval(myOpts.geom["offset"])
            spacing = "\n"*3
            print("{0}!!!WARNING: 'offset' found in BOUT.inp, "
                  "running as annulus!!!{0}".format(spacing))
            self._rho = offset + dx * np.array(np.arange(0.5, innerPoints))
        except KeyError:
            # This is the default
            self._rho = dx * np.array(np.arange(0.5, innerPoints))

        if xguards:
            # Insert the first and last grid point
            self._rho = np.insert(self._rho, 0, - 0.5*dx)
            self._rho = np.append(self._rho, self._rho[-1] + dx)
        #}}}

        #{{{z
        dy  = collect("dy"             ,\
                      path    = path   ,\
                      xguards = xguards,\
                      yguards = yguards,\
                      info    = False)
        MYG = collect("MYG"            ,\
                      path    = path   ,\
                      xguards = xguards,\
                      yguards = yguards,\
                      info    = False)

        nPoints  = dy.shape[1]
        self._dy = dy[0,0]

        if yguards:
            innerPoints = nPoints - 2*MYG
        else:
            innerPoints = nPoints

        self._z = self._dy * np.array(np.arange(0.5, innerPoints))

        if yguards:
            # Insert the first and last grid point
            self._z = np.insert(self._z, 0, - 0.5*self._dy)
            self._z = np.append(self._z, self._z[-1] + self._dy)
        #}}}

        #{{{theta
        self._dz = collect("dz"             ,\
                           path    = path   ,\
                           xguards = xguards,\
                           yguards = yguards,\
                           info    = False)
        MZ       = collect("MZ"             ,\
                           path    = path   ,\
                           xguards = xguards,\
                           yguards = yguards,\
                           info    = False)

        # Subtract the unused plane
        innerPoints = MZ - 1

        self._theta = self._dz * np.array(np.arange(0.0, innerPoints))
        #}}}

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
        self._t =\
            collect("t_array", path=self._path, tind=self._tind, info=False)

        # Slice in t
        if self._tSlice is not None:
            if self._tSlice.step is not None:
                self._t = self._t[::self._tSlice.step]

        # Convert to physical units
        if self.physicalU:
            try:
                self._omCI = collect("omCI", path=self._path, info=False)
                self._rhoS = collect("rhoS", path=self._path, info=False)
                self._n0   = collect("n0"  , path=self._path, info=False)
                self._Te0  = collect("Te0" , path=self._path, info=False)
                self._Ti0  = collect("Ti0" , path=self._path, info=False)
                self._B0   = collect("B0"  , path=self._path, info=False)
                self._Sn   = collect("Sn"  , path=self._path, info=False)

            except ValueError:
                # An OSError is thrown if the file is not found
                message = ("{0}{1}WARNING: Normalized quantities not found. "
                           "The time remains normalized".format("\n"*3,"!"*3))
                print(message)

                # Reset physicalU
                self.physicalU = False

        # Convert to physical units
        if self.physicalU:
            self._t   /= self._omCI
            self._rho *= self._rhoS
            self._z   *= self._rhoS

        self._frames = len(self._t)

        # Set subPolAvg option
        self._subPolAvg = subPolAvg

        # Prepare the labels
        if self.physicalU:
            self._timeTxt   = r"$t =${} $s$"
            self._rhoPosTxt = r"$\rho$ $[m]$"
            self._zPosTxt   = r"$z$ $[m]$"
        else:
            self._timeTxt   = "$t\\omega_{{ci}} =$ {}"
            self._rhoPosTxt = r"$\rho/\rho_s$"
            self._zPosTxt   = r"$z/\rho_s$"

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
        curIndices : list
            A list of the start and the stop values in the slice
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

# FIXME: This can be generalized for all plots
    #{{{_plotNumberFormatter
    def _plotNumberFormatter(self, val, pos):
        """
        Formatting numbers in the plot

        Parameters
        ----------
        val : float
            The value.
        pos : [None | float]
            The position (needed as input from FuncFormatter).
        """

        tickString = "${:.3g}".format(val)
        if "e+" in tickString:
            tickString = tickString.replace("e+" , r"\cdot 10^{")
            tickString = tickString.replace("e+0", r"\cdot 10^{")
            tickString += "}$"
        elif "e-" in tickString:
            tickString = tickString.replace("e-" , r"\cdot 10^{-")
            tickString = tickString.replace("e-0", r"\cdot 10^{-")
            tickString += "}$"
        else:
            tickString += "$"

        return tickString
    #}}}

# FIXME: Can this be genralized more?
    #{{{_getUnitsAndSetPhysical
    def _getUnitsAndSetPhysical(self, varName, var):
        """
        Calculates physical parameters from the normalized if
        self.physicalU is set, sets self._units.

        Parameters
        ----------
        varName : str
            Name of the variable.
        var : array
            The variable.

        Returns
        -------
        var : array
            The variable after eventual processing
        """
        if self.physicalU:
            # Calculate back to physical units
            if varName == "n":
                var *= self._n0
                self._units = r"\mathrm{m}^{-3}"
            elif varName == "vort":
                var *= self._omCI
                self._units = r"\mathrm{s}^{-1}"
            elif varName == "vortD":
                var *= self._omCI*self._n0
                self._units = r"\mathrm{m}^{-3}\mathrm{s}^{-1}"
            elif varName == "phi":
                var *= self._Te0/cst.e
                self._units = r"\mathrm{J}\mathrm{C}^{-1}"
            elif varName == "jPar":
                var *= cst.m_p*self._n0*self._rhoS*self._omCI
                self._units = r"\mathrm{C}\mathrm{s}^{-1}"
            elif varName == "momDensPar":
                # momDensPar is divided by m_i, so we need to multiply
                # by m_i again here
                var *= cst.m_p*self._rhoS*self._omCI*self._n0
                self._units = r"\mathrm{kg }\mathrm{m}^{-2}\mathrm{s}^{-1}"
            elif varName == "uIPar":
                var *= self._rhoS*self._omCI
                self._units = r"\mathrm{m}\mathrm{s}^{-1}"
            elif varName == "uEPar":
                var *= self._rhoS*self._omCI
                self._units = r"\mathrm{m}\mathrm{s}^{-1}"
            elif varName == "S":
                var *= self._omCI*self._n0
                self._units = r"\mathrm{m}^{-3}\mathrm{s}^{-1}"
            else:
                self._units = " "
        else:
            # Calculate back to physical units
            if varName == "n":
                self._units = r"/n_0"
            elif varName == "vort":
                self._units = r"/\omega_{{ci}}"
            elif varName == "vortD":
                self._units = r"/\omega_{{ci}}n_0"
            elif varName == "phi":
                self._units = r" q/T_{{e,0}}"
            elif varName == "jPar":
                self._units = r"/n_0c_sq"
            elif varName == "momDensPar":
                # by m_i again here
                self._units = r"/m_in_0c_s"
            elif varName == "uIPar":
                self._units = r"/c_s"
            elif varName == "uEPar":
                self._units = r"/c_s"
            elif varName == "S":
                self._units = r"/\omega_{{ci}}n_0"
            else:
                self._units = " "

        return var
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
        # FIXME: Check if this can go to the parent constructor
        #{{{x-direction
        if kwargs["xSlice"] == slice(0,None):
            self._xAx = self._rho

            # Set the label and the title
            if self.physicalU:
                self._xlabel = r"$\rho$ $[m]$"
                self._title  = r"$\theta={:.2g}^{{\circ}},$ $z={:.2g}$ $m$  "
            else:
                self._xlabel = r"$\rho/\rho_s$"
                self._title  = r"$\theta={:.2g}^{{\circ}},$ $z/\rho_s={:.2g}$ "

            self._title = self._title.\
                            format(self._theta[kwargs["zSlice"]]*(180/np.pi),\
                                   self._z    [kwargs["ySlice"]])
            # Set direction (used in save)
            self._direction = "radial"
        #}}}

        #{{{y-direction
        if kwargs["ySlice"] == slice(0,None):
            self._xAx = self._z

            # Set the label and the title
            if self.physicalU:
                self._xlabel = r"$z$ $[m]$"
                self._title  = r"$\rho={:.2g}$ $m$, $\theta={:.2g}^{{\circ}}$ "
            else:
                self._xlabel = r"$z/\rho_s$"
                self._title  =\
                    r"$\rho/\rho_s={:.2g}$, $\theta={:.2g}^{{\circ}}$ "

            self._title = self._title.\
                        format(self._rho  [kwargs["xSlice"]],\
                               self._theta[kwargs["zSlice"]]*(180/np.pi))

            # Set direction (used in save)
            self._direction = "parallel"
        #}}}

        #{{{z-direction
        if kwargs["zSlice"] == slice(0,None):
            self._xAx = self._theta

            # Set the label and the title
            if self.physicalU:
                self._xlabel = r"$\theta$"
                self._title  = r"$\rho={:.2g}$ $m$, $z={:.2g}$ $m$  "
            else:
                self._xlabel = r"$\theta$"
                self._title  = r"$\rho/\rho_s={:.2g}$, $z/rho_s={:.2g}$  "

            self._title = self._title.\
                        format(self._rho[kwargs["xSlice"]],\
                               self._z  [kwargs["ySlice"]])

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
        timeString = self._plotNumberFormatter(self._t[tInd], None)
        curTimeTxt = self._timeTxt.format(timeString)
        fig.suptitle(self._title + curTimeTxt)
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
                FuncFormatter(self._plotNumberFormatter)\
                                                   )
# FIXME: Generalize
            line.ax.get_xaxis().set_major_formatter(\
                FuncFormatter(lambda val, pos:"${:.2g}$".format(val))\
                                                   )
            # Set grid
            line.ax.grid(b = True)

        if orgObj.useCombinedPlot:
            # Find the max and the min
            allMax = np.max(allMax)
            allMin = np.min(allMin)

            # Set the y-axis limits
            orgObj.combLine.ax.set_ylim(allMin, allMax)

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
                pass

            # If Variable not saved each timestep
            if len(line.field.shape) == 3:
                # Make it a 4d variable
                field      = np.zeros(( len(self._t), *line.field.shape))
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
                pass

            # If Variable not saved each timestep
            if len(line.field.shape) == 3:
                # Make it a 4d variable
                field      = np.zeros(( len(self._t), *line.field.shape))
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
            self.physicalU = False

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
                # Save the animation
                anim.save(fileName + ".gif"              ,\
                          writer         = "imagemagick"   ,\
                          savefig_kwargs = {"pad_inches":0},\
                          )
                print("Saved to {}.gif".format(fileName))
        else:
            if self._savePlot:
                # Save the figure
                self._fig.savefig(fileName + ".pdf"  ,\
                                  transparent = False  ,\
                                  bbox_inches = "tight",\
                                  pad_inches  = 0      ,\
                                  )
                print("Saved to {}.pdf".format(fileName))

        if self._showPlot:
            self._fig.show()

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
    def __init__(self                      ,\
                 path                      ,\
                 varName                   ,\
                 var        = None         ,\
                 xguards    = False        ,\
                 yguards    = False        ,\
                 varMax     = None         ,\
                 varMin     = None         ,\
                 varyMaxMin = False        ,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Builds the mesh from the CylinderMesh class
        2. Either:
            1. Collects the variable
            2. Calculates the variable form varFunc if set
        3. Calculates the physical units if physicalU is set
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
        varMax : [None | float]
            Setting a hard upper limit z-axis in the plot.
        varMin : [None | float]
            Setting a hard lower limit z-axis in the plot.
        varyMaxMin : bool
            Whether or not the limits of the z-axis should be
            set to the max/min of the current timestep or not.
        xguards : bool
            If xguards was used when collecting (used to set the mesh).
        yguards : bool
            If yguards was used when collecting (used to set the mesh).
        **kwargs : keyword arguments
            See the constructor of Plot for details.
        """
        #}}}

        # Call the constructor of the parent class
        super(Plot2D, self).__init__(path   ,\
                                     xguards,\
                                     yguards,\
                                     **kwargs)

        # Check that the indices are properly set
        if (kwargs["xSlice"] == slice(0,None)) and\
           (kwargs["ySlice"] == slice(0,None)) and\
           (kwargs["zSlice"] == slice(0,None)):
            message = "3 slices were given, although only 2 is possible"
            raise ValueError(message)

        # Make it possible to filter warnings (f.ex if no variation in the data)
        warnings.filterwarnings("error")

        # Set member data from the index
        self._varyMaxMin = varyMaxMin
        self._varMax     = varMax
        self._varMin     = varMin

        # Set additional plot properties
        self._latexSize = 35
        self._nCont     = 100
        self._pltName   = None

        # Create a CylinderMesh object
        self._cyl = CylinderMesh(self._rho, self._theta, self._z, xguards)

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
                self._cyl.addLastThetaSlice(self._variable, len(self._t))

        self._variable =\
                self._getUnitsAndSetPhysical(varName, self._variable)

        if xguards:
            # Remove the inner ghost points from the variable
            self._variable = np.delete(self._variable, (0), axis=1)

        # Get the max and the min so that we can keep the color coding correct
        if self._varMax == None:
            self._varMax = np.max(self._variable)
        if self._varMin == None:
            self._varMin = np.min(self._variable)

        # We need to manually sepcify the levels in order to have a
        # fixed color bar
        self._levels = np.linspace(self._varMin  ,\
                                   self._varMax  ,\
                                   self._nCont   ,\
                                   endpoint = True)

        # Then theta index corresponding to pi
        piInd = round(self._variable.shape[3]/2)

        # Get the Z values of the X, Y, Z plots
        # We subscript the last index of self._@ind, as this is given as
        # a range in the Plot constructor
        self._Z_RT = self._variable[:, :, self._yind[-1], :             ]
        self._Z_RZ = self._variable[:, :, :             , self._zind[-1]]
        # Get the Z value in the RZ plane which is pi from the current index
        if self._zind[-1] > piInd:
            self._Z_RZ_P_PI = self._variable[:, :, :, self._zind[-1] - piInd]
        else:
            self._Z_RZ_P_PI = self._variable[:, :, :, self._zind[-1] + piInd]

        self._setLines()

        # Create the figure and axis
        pltSize      = (30,15)
        gs           = GridSpec(1, 3, width_ratios=[20, 20, 1])
        self._fig    = plt.figure(figsize = pltSize)
        self._ax1    = self._fig.add_subplot(gs[0])
        self._ax2    = self._fig.add_subplot(gs[1])
        self._cBarAx = self._fig.add_subplot(gs[2])
        self._fig.subplots_adjust(wspace=0.25)
        self._ax1.grid(True)
        self._ax2.grid(True)

        # Create placeholder for colorbar and images
        self._cbarPlane = None
        self._images = []
        #}}}

    #{{{setLines
    def _setLines(self):
        """ Set the lines which shows where the data is sliced"""

        # The slice lines we are plotting
        rhoStart = self._rho[0]
        rhoEnd   = self._rho[-1]

        # Calculate the numerical value of the theta angle and the z value
        thetaRad         = self._dz*self._zind[-1]
        thetaPPi         = thetaRad + np.pi
        self._zVal       = self._z[self._ySlice]
        self._thetaDeg   = thetaRad*(180/np.pi)

        # Set coordinates for the lines which indicate how the data is
        # sliced
        # Organized in start and stop pairs
        # We need two lines due to the center of the cylinder
        self._RTLine1XVals=(rhoStart*np.cos(thetaRad), rhoEnd*np.cos(thetaRad))
        self._RTLine1YVals=(rhoStart*np.sin(thetaRad), rhoEnd*np.sin(thetaRad))
        self._RTLine2XVals=(rhoStart*np.cos(thetaPPi), rhoEnd*np.cos(thetaPPi))
        self._RTLine2YVals=(rhoStart*np.sin(thetaPPi), rhoEnd*np.sin(thetaPPi))
        self._RZLine1XVals=(-rhoEnd                  , -rhoStart              )
        self._RZLine1YVals=(self._zVal               , self._zVal             )
        self._RZLine2XVals=(rhoStart                 , rhoEnd                 )
        self._RZLine2YVals=(self._zVal               , self._zVal             )
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

        # If we want the max and min to vary
        if self._varyMaxMin:
            # Update the max and min
            self._varMax = np.max(self.Z_RT)
            self._varMin = np.min(self.Z_RT)
            # Update the levels
            self._levels = np.linspace(self._varMin   ,\
                                       self._varMax   ,\
                                       self._nCont    ,\
                                       endpoint = True,\
                                       )

        # Check that levels are rising
        if not(self._levels is None):
            if len(self._levels) > 1 and np.amin(np.diff(self._levels)) <= 0.0:
                self._levels = None

        Z_RT      = self._Z_RT     [tInd, :, :]
        Z_RZ      = self._Z_RZ     [tInd, :, :]
        Z_RZ_P_PI = self._Z_RZ_P_PI[tInd, :, :]

        # Plot the perpendicular plane
        perpPlane  = self._ax1.contourf(self._cyl.X_RT       ,\
                                        self._cyl.Y_RT       ,\
                                        Z_RT                 ,\
                                        cmap   = cm.RdYlBu_r ,\
                                        vmax   = self._varMax,\
                                        vmin   = self._varMin,\
                                        levels = self._levels,\
                                        )

        # Plot the parallel plane
        parPlane  = self._ax2.contourf(self._cyl.X_RZ       ,\
                                       self._cyl.Y_RZ       ,\
                                       Z_RZ                 ,\
                                       cmap   = cm.RdYlBu_r ,\
                                       vmax   = self._varMax,\
                                       vmin   = self._varMin,\
                                       levels = self._levels,\
                                       )
        parPlaneNeg  = self._ax2.contourf(self._cyl.X_RZ_NEG   ,\
                                          self._cyl.Y_RZ       ,\
                                          Z_RZ_P_PI            ,\
                                          cmap   = cm.RdYlBu_r ,\
                                          vmax   = self._varMax,\
                                          vmin   = self._varMin,\
                                          levels = self._levels,\
                                          )

        # Draw the grids
        self._ax1.grid(b=True)
        self._ax2.grid(b=True)

        # x and y labels
        self._ax1.set_xlabel(self._rhoPosTxt, fontsize = self._latexSize)
        self._ax1.set_ylabel(self._rhoPosTxt, fontsize = self._latexSize)
        self._ax2.set_xlabel(self._rhoPosTxt, fontsize = self._latexSize)
        self._ax2.set_ylabel(self._zPosTxt  , fontsize = self._latexSize)


        # Title preparation
        if self.physicalU:
            fixedParTxt = "{:.2g}$ $m".format(self._z[self._ySlice])
        else:
            fixedParTxt = "{:.2g}".format(self._z[self._ySlice])
        timeString = self._plotNumberFormatter(self._t[tInd], None)
        timeTxt = self._timeTxt.format(timeString)

        # Titles
        ax1Title = r"$z= {}$".format(fixedParTxt)
        ax2Title = r"$\theta={:.0f}^{{\circ}}$".format(self._thetaDeg)

        # Title axis 1
        self._ax1txt = self._ax1.text(0.5, 1.05,\
                                      ax1Title,\
                                      horizontalalignment = "center",\
                                      verticalalignment = "center",\
                                      fontsize = self._latexSize,\
                                      transform = self._ax1.transAxes)

        # Title axis 2
        self._ax2txt = self._ax2.text(0.5, 1.05,\
                            ax2Title,\
                            horizontalalignment = "center",\
                            verticalalignment = "center",\
                            fontsize = self._latexSize,\
                            transform = self._ax2.transAxes)

        # Title mid
        # Text for the figure. Could append this to the figure itself,
        # but it seems to be easier to just add it to an axis due to
        # animation
        self._figTxt = self._ax1.text(1.10, 1.05,\
                                      timeTxt,\
                                      horizontalalignment = "center",\
                                      verticalalignment = "center",\
                                      fontsize = self._latexSize,\
                                      transform = self._ax1.transAxes)



        self._ax1.get_xaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:"${:.2g}$".format(val))\
                                                   )
        self._ax1.get_yaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:"${:.2g}$".format(val))\
                                                   )
        self._ax2.get_xaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:"${:.2g}$".format(val))\
                                                   )
        self._ax2.get_yaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:"${:.2g}$".format(val))\
                                                   )

        # Make the axis equal
        # FIXME: Make it an option to not plot ax2 equal => can see ||
        #        structures
        self._ax1.axis("equal")
        self._ax2.axis("equal")

        # Current API inconsistency fix
        # (https://github.com/matplotlib/matplotlib/issues/6139):
        addArtPerpPlane   = perpPlane.collections
        addArtParPlane    = parPlane.collections
        addArtParPlaneNeg = parPlaneNeg.collections

        # Put images together
        self._images.append(addArtPerpPlane + [self._ax1txt] +\
                            addArtParPlane + addArtParPlaneNeg +\
                            [self._ax2txt] + [self._figTxt])

        if self._cbarPlane is None:
            self._cbarPlane = parPlane
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

        # The colorbar needs only to be plotted once
        # Make the colorbar
        # format = "%.g" gave undesired results
        try:
            cbar = self._fig.colorbar(self._cbarPlane                   ,\
                                      cax    = self._cBarAx             ,\
                                      format = FuncFormatter(     \
                                              self._plotNumberFormatter),\
                                      )
            if self.physicalU:
                cbarName = r"${}$ $[{}]$".format(self._pltName, self._units)
            else:
                cbarName = r"${}{}$".format(self._pltName, self._units)

            cbar.set_label(label = cbarName, size = titleSize + 5)

        except RuntimeWarning:
            message  = "RuntimeError caught in cbar in " + self._pltName
            message += ". No cbar will be set!"
        # Lines needs only to be plotted once
        # Par line 1
        self._ax2.plot(self._RZLine1XVals,\
                       self._RZLine1YVals,\
                       "--k"             ,\
                       linewidth = 1     ,\
                       )
        # Par line 2
        self._ax2.plot(self._RZLine2XVals,\
                       self._RZLine2YVals,\
                       "--k"             ,\
                       linewidth = 1     ,\
                       )
        # Perp line 1
        self._ax1.plot(self._RTLine1XVals,\
                       self._RTLine1YVals,\
                       "--k"             ,\
                       linewidth = 1     ,\
                       )
        # Perp line 2
        self._ax1.plot(self._RTLine2XVals,\
                       self._RTLine2YVals,\
                       "--k"             ,\
                       linewidth = 1     ,\
                       )

        # Need to specify rect in order to have top text
        self._fig.tight_layout(w_pad = 2.5, rect=[0,0,1,0.97])

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
            fileName = saveName + "-2D"
            fileName = os.path.join(savePath, fileName)

        # Animate if we have more than one frame
        if self._frames > 1:
            # Create the plots
            for tInd in range(1, self._frames):
                self._plot2D(tInd)

            # Animate
            anim = animation.ArtistAnimation(self._fig            ,\
                                             self._images         ,\
                                             blit   = False       ,\
                                             )

            if self._savePlot:
                # Save the animation
                anim.save(fileName + ".gif"              ,\
                          writer         = "imagemagick"   ,\
                          savefig_kwargs =\
                            {"pad_inches"         :0       ,\
                             "bbox_extra_artists" :(cbar,\
                                                   self._ax1txt,\
                                                   self._ax2txt),\
                             "bbox_inches"        :"tight"},\
                          )
                print("Saved to {}.gif".format(fileName))
        else:
            if self._savePlot:
                # Save the figure
                self._fig.savefig(fileName + ".pdf"  ,\
                                  transparent = False  ,\
                                  bbox_inches = "tight",\
                                  pad_inches  = 0      ,\
                                  )
                print("Saved to {}.pdf".format(fileName))

        if self._showPlot:
            self._fig.show()

        plt.close(self._fig)
    #}}}
#}}}
