#!/usr/bin/env python

"""
Contains class for plotting
"""

from ..statistics import polAvg
from .getStrings import getSaveString
from .cylinderMesh import CylinderMesh
from matplotlib import get_backend
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from boutdata import collect
from boututils.options import BOUTOptions
import numpy as np
import warnings

# All post processing functions called by bout_runners must accept the
# first argument from bout_runners (called 'folder' in
# __call_post_processing_function)

#{{{class Plot
class Plot(object):
    """
    Parent class for plotting the results of the CELMA code.

    Handles:

    * Setting the collect options
    * Formation of ticks
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
                 polAvg     = False        ,\
                 showPlot   = False        ,\
                 savePlot   = True         ,\
                 saveFolder = None         ,\
                ):
        #{{{docstring
        """
        The constructor sets the member data

        Input:
        path         - The path to collect from
        xguards      - If xguards should be included when collecting
        yguards      - If yguards should be included when collecting
        xSlice       - How the data will be sliced in x
        ySlice       - How the data will be sliced in y
        zSlice       - How the data will be sliced in z
        tSlice       - How the data will be sliced in t
        polAvg       - Whether or not to perform a poloidal average of
                       the data
        showPlot     - If the plot should be displayed
        savePlot     - If plot should be saved
        saveFolder   - Name of the folder to save plots in
        """
        #}}}

        # Set member data from input
        self._path       = path
        self._xguards    = xguards
        self._yguards    = yguards
        self._showPlot   = showPlot
        self._savePlot   = savePlot
        self._saveFolder = saveFolder

        #{{{ Set the plot style
        self._titleSize = 30
        plt.rc("font",   size      = 30)
        plt.rc("axes",   labelsize = 25, titlesize = self._titleSize)
        plt.rc("xtick",  labelsize = 25)
        plt.rc("ytick",  labelsize = 25)
        plt.rc("legend", fontsize  = 20)
        plt.rc("lines",  linewidth = 2)
        #}}}

        # Get the coordinates
        #{{{rho
        dx = collect('dx'             ,\
                     path    = path   ,\
                     xguards = xguards,\
                     yguards = yguards,\
                     info    = False)
        MXG = collect('MXG'            ,\
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
            offset = eval(myOpts.geom['offset'])
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
        dy  = collect('dy'             ,\
                      path    = path   ,\
                      xguards = xguards,\
                      yguards = yguards,\
                      info    = False)
        MYG = collect('MYG'            ,\
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
        self._dz = collect('dz'             ,\
                           path    = path   ,\
                           xguards = xguards,\
                           yguards = yguards,\
                           info    = False)
        MZ       = collect('MZ'             ,\
                           path    = path   ,\
                           xguards = xguards,\
                           yguards = yguards,\
                           info    = False)

        # Subtract the unused plane
        innerPoints = MZ - 1

        self._theta = self._dz * np.array(np.arange(0.0, innerPoints))
        #}}}

        # Get proper indices
        self._xind = self._getIndices(xSlice)
        self._yind = self._getIndices(ySlice)
        self._zind = self._getIndices(zSlice)
        self._tind = self._getIndices(tSlice)
        # Used if we are taking poloidal averages
        self._xSlice = xSlice
        self._ySlice = ySlice
        self._zSlice = zSlice

        # Get the time
        self._t = collect('t_array', path=self._path, tind=self._tind, info=False)

        self._frames = len(self._t)

        # Set polAvg option
        self._polAvg = polAvg
    #}}}

    #{{{ _getIndices
    def _getIndices(self, curSlice):
        """
        Return the slice such that it can be given as an input to 'collect'
        """

        if type(curSlice) == slice:
            curIndices = []
            curIndices.append(curSlice.start)
            if curSlice.stop == None:
                curIndices = None
            else:
                curIndices.append(curSlice.stop)
        elif curSlice is None:
            curIndices = curSlice
        else:
            curIndices = [curSlice, curSlice]

        return curIndices
    #}}}

    #{{{_plotNumberFormatter
    def _plotNumberFormatter(self, val, pos):
        """
        Formatting numbers in the plot

        Input
        val - The value
        pos - The position (needed as input from FuncFormatter)
        """

        tickString = '${:g}'.format(val)
        tickString = tickString.replace('e+', r'\cdot 10^{')
        tickString = tickString.replace('e-', r'\cdot 10^{-')
        tickString += '}$'

        return tickString
    #}}}
#}}}

#{{{class Plot1D
class Plot1D(Plot):
    """
    Class for plotting the results of the CELMA code in 1D.
    Inherits from the Plot class

    Handles:

    * Collection of the variables
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
        The constructor sets the member data

        Input specific for Plot1D:
        marker  - The type of marker to be used in the plot

        For the other input, refer to the docstring of the Plot class
        """
        #}}}

        # Call the constructor of the parent class
        super(Plot1D, self).__init__(*args, **kwargs)

        # Check that the indices are properly set
        # Note that this is set after super, as super will check for bad
        # input
        if (kwargs['xSlice'] == slice(0,None)) and\
           (kwargs['ySlice'] == slice(0,None)) and\
           (kwargs['zSlice'] == slice(0,None)):
            message = "3 slices were given, although only 1 is possible"
            raise RuntimeError(message)
        elif (kwargs['xSlice'] == slice(0,None) and\
              kwargs['ySlice'] == slice(0,None)) or\
             (kwargs['ySlice'] == slice(0,None) and\
              kwargs['zSlice'] == slice(0,None)) or\
             (kwargs['zSlice'] == slice(0,None) and\
              kwargs['xSlice'] == slice(0,None)):
            message = "2 slices were given, although only 1 is possible"
            raise RuntimeError(message)

        # Get the x-axis of the plot
        #{{{x-direction
        if kwargs['xSlice'] == slice(0,None):
            self._xAx = self._rho

            # Set the label and the title
            self._xlabel = r'$\rho$'
            self._title  = r'$\theta_i={0},$ $z_i={1}$  '.\
                    format(kwargs['zSlice'], kwargs['ySlice'])

            # Set direction (used in save)
            self._direction = 'radial'
        #}}}

        #{{{y-direction
        if kwargs['ySlice'] == slice(0,None):
            self._xAx = self._z

            # Set the label and the title
            self._xlabel = r'$z$'
            self._title  = r'$\rho_i={0}$, $\theta_i={1}$  '.\
                    format(kwargs['xSlice'], kwargs['zSlice'])

            # Set direction (used in save)
            self._direction = 'parallel'
        #}}}

        #{{{z-direction
        if kwargs['zSlice'] == slice(0,None):
            self._xAx = self._theta

            # Set the label and the title
            self._xlabel = r'$\theta$'
            self._title  = r'$\rho={0}$, $z_i={1}$  '.\
                    format(kwargs['xSlice'], kwargs['ySlice'])

            # Set direction (used in save)
            self._direction = 'theta'
        #}}}

        # Set the input data
        self._marker = marker
    #}}}

    #{{{_animFunction
    def _animFunction(self, tInd, orgObj, fig):
        """
        Function which updates the data.

        As blitting is False, there is no reason to return the lines
        """

        # Plot the lines
        for ind, line in enumerate(orgObj.lines):
            # Plot the normal lines
            line.lineObj.set_data(self._xAx, line.field[tInd,:])

            if orgObj.useCombinedPlot:
                # Plot the line in the combined plot
                orgObj.combLineLineObjs[ind].\
                        set_data(self._xAx, line.field[tInd,:])

        timeString = self._plotNumberFormatter(self._t[tInd], None)
        fig.suptitle(self._title + r'$\omega_{ci}^{-1} = $' + timeString)
    #}}}

    #{{{_plotLines
    def _plotLines(self, fig, orgObj, tInd):
        """
        Plots the other lines into the combined line plot

        Input
        fig    - The figure
        orgObj - Organizer object
        tInd   - The time index to plot for
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
                line.ax.tick_params(labelbottom='off')

            line.ax.legend(loc='upper right', fancybox=True, framealpha=0.5, numpoints=1)
            # Avoid ticks collision
            line.ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
            line.ax.locator_params(axis='y', tight=True, nbins=6)
            # Avoid silly top value
            line.ax.get_yaxis().get_major_formatter().set_useOffset(False)
            # Use own fuction to deal with ticks
            line.ax.get_yaxis().set_major_formatter(\
                FuncFormatter(self._plotNumberFormatter)\
                                                   )
            line.ax.get_xaxis().set_major_formatter(\
                FuncFormatter(lambda val, pos:'${:d}$'.format(int(val)))\
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
        if get_backend() == 'QT4Agg':
            figManager = plt.get_current_fig_manager()
            figManager.window.showMaximized()
    #}}}

    #{{{collectLine
    def collectLine(self, line):
        """Collects the data for one line and reshapes it"""

        if self._polAvg:
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
                # FIXME: What is the point of this?
                # I just added the below
                line.field = collect(line.name,\
                                     path    = self._path   ,\
                                     xguards = self._xguards,\
                                     yguards = self._yguards,\
                                     tind    = self._tind   ,\
                                     info    = False)
                # pass

            # If Variable not saved each timestep
            if len(line.field.shape) == 3:
                # Make it a 4d variable
                field      = np.zeros(( len(self._t), *line.field.shape))
                # Copy the field in to each time
                field[:]   = line.field
                line.field = field

            # Take the poloidal average, and slice the result
            line.field = polAvg(line.field) \
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

        # Flatten the variables except the time dimension
        # -1 => total size divided by product of all other listed dimensions
        line.field = line.field.reshape(line.field.shape[0], -1)
    #}}}

    #{{{plotDriver
    def plotDriver(self, fig, orgObj, timeFolder):
        """
        Function which drives the plotting.

        Input
        fig        - The figure
        orgObj     - The organization object
        timeFolder - Name of the timeFolder (if none is given, one is
                     going to be made)

        Output
        timeFolder - The timefolder used when eventually saving the plot
        """

        # Initial plot
        self._plotLines(fig, orgObj, 0)

        if self._savePlot:
            # Make a saveName by stripping the orgObj's plot name for bad
            # characters
            saveName = orgObj.pltName.replace("\\", "")
            saveName = saveName.replace("{", "")
            saveName = saveName.replace("}", "")
            saveName = saveName.replace("^", "")
            fileName = saveName + '-' + self._direction
            prePaths = ['visualization', self._saveFolder]
            if self._polAvg:
                postPaths = 'polAvg'
            else:
                postPaths = []
            saveString, timeFolder = getSaveString(fileName               ,\
                                                   self._path             ,\
                                                   timeFolder = timeFolder,\
                                                   prePaths   = prePaths  ,\
                                                   postPaths  = postPaths ,\
                                                   )

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
                anim.save(saveString + '.gif'              ,\
                          writer         = 'imagemagick'   ,\
                          savefig_kwargs = {'pad_inches':0},\
                          )
                print("Saved to {}.gif".format(saveString))
        else:
            if self._savePlot:
                # Save the figure
                plt.savefig(saveString + '.png'  ,\
                            transparent = False  ,\
                            bbox_inches = 'tight',\
                            pad_inches  = 0      ,\
                            )
                print("Saved to {}.png".format(saveString))

        if self._showPlot:
            plt.show()

        return timeFolder
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
                 xguards    = False        ,\
                 yguards    = False        ,\
                 xSlice     = slice(0,None),\
                 ySlice     = slice(0,None),\
                 zSlice     = slice(0,None),\
                 varMax     = None         ,\
                 varMin     = None         ,\
                 varyMaxMin = False        ,\
                 varFunc    = None         ,\
                 **kwargs):
        #{{{docstring
        """
        The constructor sets the member data

        Specific input for Plot2D
        varName    - Name of the field which is going to be collected
        varMax     - Setting a hard upper limit z-axis in the plot
        varMin     - Setting a hard lower limit z-axis in the plot
        varyMaxMin - Whether or not the limits of the z-axis should be
                     set to the max/min of the current timestep or not
        xguards    - If xguards should be included when collecting
        yguards    - If yguards should be included when collecting
        xSlice     - How the data will be sliced in x
        ySlice     - How the data will be sliced in y
        zSlice     - How the data will be sliced in z
        varFunc    - Function which returns the variable (used if
                     variables is not collectable)
        For more details, refer to the docstring of the Plot class
        """
        #}}}

        # Call the constructor of the parent class
        super(Plot2D, self).__init__(path   ,\
                                     xguards,\
                                     yguards,\
                                     xSlice ,\
                                     ySlice ,\
                                     zSlice ,\
                                     **kwargs)

        # Check that the indices are properly set
        if (xSlice == slice(0,None)) and\
           (ySlice == slice(0,None)) and\
           (zSlice == slice(0,None)):
            message = "3 slices were given, although only 2 is possible"
            raise RuntimeError(message)

        # Make it possible to filter warnings (f.ex if no variation in the data)
        warnings.filterwarnings('error')

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
        if varFunc is None:
            self._variable = collect(varName             ,\
                                     path    = path      ,\
                                     yguards = yguards   ,\
                                     xguards = xguards   ,\
                                     tind    = self._tind,\
                                     info    = False     ,\
                                     )
        else:
            self._variable = varFunc(path    = path      ,\
                                     yguards = yguards   ,\
                                     xguards = xguards   ,\
                                     tind    = self._tind,\
                                     info    = False     ,\
                                     **kwargs)

        if self._polAvg:
            self._variable = polAvg(self._variable)

        # Add the last theta slice
        self._variable =\
                self._cyl.addLastThetaSlice(self._variable, len(self._t))

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

        # The slice lines we are plotting
        rhoStart = self._rho[0]
        rhoEnd   = self._rho[-1]
        zStart   = self._z  [0]

        # Calculate the numerical value of the theta angle and the z value
        thetaRad         = self._dz*self._zind[-1]
        thetaPPi         = thetaRad + np.pi
        self._zVal       = zStart + self._yind[-1]*self._dy
        self._thetaDeg   = thetaRad*(180/np.pi)

        # Set coordinates for the lines which indicate how the data is
        # sliced
        # Organized in start and stop pairs
        # We need two lines due to the center of the cylinder
        self._RTLine1XVals = (rhoStart*np.cos(thetaRad), rhoEnd*np.cos(thetaRad))
        self._RTLine1YVals = (rhoStart*np.sin(thetaRad), rhoEnd*np.sin(thetaRad))
        self._RTLine2XVals = (rhoStart*np.cos(thetaPPi), rhoEnd*np.cos(thetaPPi))
        self._RTLine2YVals = (rhoStart*np.sin(thetaPPi), rhoEnd*np.sin(thetaPPi))
        self._RZLine1XVals = (-rhoEnd                  , -rhoStart              )
        self._RZLine1YVals = (self._zVal               , self._zVal             )
        self._RZLine2XVals = (rhoStart                 , rhoEnd                 )
        self._RZLine2YVals = (self._zVal               , self._zVal             )

        # Create the figure and axis
        pltSize      = (30,15)
        gs           = GridSpec(1, 3, width_ratios=[20, 20, 1])
        self._fig    = plt.figure(figsize = pltSize)
        self._ax1    = plt.subplot(gs[0])
        self._ax2    = plt.subplot(gs[1])
        self._cBarAx = plt.subplot(gs[2])
        self._fig.subplots_adjust(wspace=0.25)
        self._ax1.grid(True)
        self._ax2.grid(True)
        #}}}

    #{{{_plot2D
    def _plot2D(self, tInd):
        #{{{docstring
        """
        Performs the actual plotting

        Input:
        tInd - The index to plot for
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
        perpLine1 = self._ax1.plot(self._RTLine1XVals,\
                                   self._RTLine1YVals,\
                                   '--k'             ,\
                                   linewidth = 1     ,\
                                   )
        perpLine2 = self._ax1.plot(self._RTLine2XVals,\
                                   self._RTLine2YVals,\
                                   '--k'             ,\
                                   linewidth = 1     ,\
                                   )

        # Draw the grids
        self._ax1.grid(b=True)

        # Decorations
        self._ax1.set_xlabel(r'$\rho$', fontsize = self._latexSize)
        self._ax1.set_ylabel(r'$\rho$', fontsize = self._latexSize)
        self._ax1.set_title (r'$\omega_{ci}^{-1} =' + '{:g}'.format(self._t[tInd]) +\
                             r' \quad z=' +\
                             '{:.2f}'.format(self._zVal) + r'$',\
                             fontsize = self._latexSize)

        # Plot the parallel plane
        parPlane  = self._ax2.contourf(self._cyl.X_RZ       ,\
                                       self._cyl.Y_RZ       ,\
                                       Z_RZ                 ,\
                                       cmap   = cm.RdYlBu_r ,\
                                       vmax   = self._varMax,\
                                       vmin   = self._varMin,\
                                       levels = self._levels,\
                                       )
        parPlane  = self._ax2.contourf(self._cyl.X_RZ_NEG   ,\
                                       self._cyl.Y_RZ       ,\
                                       Z_RZ_P_PI            ,\
                                       cmap   = cm.RdYlBu_r ,\
                                       vmax   = self._varMax,\
                                       vmin   = self._varMin,\
                                       levels = self._levels,\
                                       )
        parLine1 = self._ax2.plot(self._RZLine1XVals,\
                                  self._RZLine1YVals,\
                                  '--k'             ,\
                                  linewidth = 1     ,\
                                  )
        parLine2 = self._ax2.plot(self._RZLine2XVals,\
                                  self._RZLine2YVals,\
                                  '--k'             ,\
                                  linewidth = 1     ,\
                                  )

        # Draw the grids
        self._ax2.grid(b=True)

        # Decorations
        self._ax2.set_xlabel(r'$\rho$', fontsize = self._latexSize)
        self._ax2.set_ylabel(r'$z$'   , fontsize = self._latexSize)
        self._ax2.set_title(r'$\omega_{ci}^{-1} =' + '{:g}'.format(self._t[tInd]) +\
                            r' \quad \theta=' +\
                           '{:.0f}'.format(self._thetaDeg) + r'^{\circ}$',\
                           fontsize = self._latexSize)

        self._ax1.get_xaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:'${:d}$'.format(int(val)))\
                                                   )
        self._ax1.get_yaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:'${:d}$'.format(int(val)))\
                                                   )
        self._ax2.get_xaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:'${:d}$'.format(int(val)))\
                                                   )
        self._ax2.get_yaxis().set_major_formatter(\
            FuncFormatter(lambda val, pos:'${:d}$'.format(int(val)))\
                                                   )

        # Make the axis equal
        self._ax1.axis('equal')
        self._ax2.axis('equal')

        # Make the colorbar
        # format = '%.g' gave undesired results
        try:
            cbar = self._fig.colorbar(parPlane                          ,\
                                      cax    = self._cBarAx             ,\
                                      format = FuncFormatter(     \
                                              self._plotNumberFormatter),\
                                      )
            cbar.set_label(r'$' + self._pltName + '$')
        except RuntimeWarning:
            message  = 'RuntimeError caught in cbar in ' + self._pltName
            message += '. No cbar will be set!'
    #}}}

    #{{{plotDriver
    def plotDriver(self, pltName, timeFolder):
        """
        Function which drived the plotting.

        Input
        pltName    - Name of the plot written in LaTeX format, but
                     without the $
        timeFolder - Name of the timeFolder (if none is given, one is
                     going to be made)

        Output
        timeFolder - The timefolder used when eventually saving the plot
        """

        self._pltName = pltName

        # Initial plot
        self._plot2D(0)

        if self._savePlot:
            # Make a saveName by stripping the orgObj's plot name for bad
            # characters
            saveName = pltName.replace("\\", "")
            saveName = saveName.replace("{", "")
            saveName = saveName.replace("}", "")
            saveName = saveName.replace("^", "")
            fileName = saveName + '-2D'
            prePaths = ['visualization', self._saveFolder]
            if self._polAvg:
                postPaths = 'polAvg'
            else:
                postPaths = []
            saveString, timeFolder = getSaveString(fileName               ,\
                                                   self._path             ,\
                                                   timeFolder = timeFolder,\
                                                   prePaths   = prePaths  ,\
                                                   postPaths  = postPaths ,\
                                                   )

        # Animate if we have more than one frame
        if self._frames > 1:
            anim = animation.FuncAnimation(self._fig            ,\
                                           self._plot2D         ,\
                                           frames = self._frames,\
                                           blit   = False       ,\
                                           )

            if self._savePlot:
                # Save the animation
                anim.save(saveString + '.gif'              ,\
                          writer         = 'imagemagick'   ,\
                          savefig_kwargs = {'pad_inches':0},\
                          )
                print("Saved to {}.gif".format(saveString))
        else:
            if self._savePlot:
                # Save the figure
                plt.savefig(saveString + '.png'  ,\
                            transparent = False  ,\
                            bbox_inches = 'tight',\
                            pad_inches  = 0      ,\
                            )
                print("Saved to {}.png".format(saveString))

        if self._showPlot:
            plt.show()

        return timeFolder
    #}}}
#}}}
