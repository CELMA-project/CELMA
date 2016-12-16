#!/usr/bin/env python

"""
Contains functions for plotting the 2D fields
"""


if self._mode == "perpAndPar".lower()\
   or self._mode == "perpAndPol".lower():
    self._createLines()

# FIXME: Remove "n" support, rather make a map!!!


# FIXME: This belongs to plotting
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

# FIXME: Put this to the plotting section
    # Set additional plot properties
    self._latexSize = 35
    self._nCont     = 100
    self._pltName   = None

    # Make it possible to filter warnings (Ex: no variation in the data)
    warnings.filterwarnings("error")


    # Get the max and the min so that we can keep the color coding correct
    if self._varMax == None:
        self._varMax = np.max(self._variable)
    if self._varMin == None:
        self._varMin = np.min(self._variable)
    # Diverging colormap for fluctuations
    if self._fluctuation:
        self._varMax = np.max([np.abs(self._varMax), np.abs(self._varMin)])
        self._varMin = - self._varMax

    # We need to manually sepcify the levels in order to have a
    # fixed color bar
    self._levels = np.linspace(self._varMin  ,\
                               self._varMax  ,\
                               self._nCont   ,\
                               endpoint = True)

    """
    Class for plotting the results of the CELMA code in 2D.
    Inherits from the Plot class

    Handles:

    * Collection of the variables
    * Plotting of the variables
    * Animation of the variables
    """

# Create the figure and axis
if self._mode == "perpAndPar".lower() or\
   self._mode == "perpAndPol".lower():
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





# FIXME: Provide max and min from the outside
#        This can be a tuple for each instance
#        Defaults to None
# FIXME: Make a getMaXMinAnimation
#        Make a getLevels2DAnimation

        # Diverging colormap for fluctuations
        if self._fluctuation:
            self._varMax =\
                    np.max([np.abs(self._varMax), np.abs(self._varMin)])
            self._varMin = - self._varMax

    # Check that levels are rising
    if not(self._levels is None):
        if len(self._levels) > 1 and np.amin(np.diff(self._levels)) <= 0.0:
            self._levels = None

# FIXME: Needed
self._levels = np.linspace(self._varMin  ,\
                           self._varMax  ,\
                           self._nCont   ,\
                           endpoint = True)

        # Update the levels
        levels = np.linspace(self._varMin   ,\
                             self._varMax   ,\
                             self._nCont    ,\
                             endpoint = True,\
                             )
    # If we want the max and min to vary
    if self._varyMaxMin and tInd:
        # Update the max and min
        self._varMax = np.max(maxList)
        self._varMin = np.min(minList)

        # Update the levels just if there is any difference
        if np.amin(np.diff(levels)) > 0.0:
            self._levels = levels




# FIXME: Make these functions: In ax, out ax
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

# FIXME: Put this on top of decorator
    # FIXME: Slice in the variable
    Z_ZT = self._Z_ZT[tInd, :, :]


    # Specify repeated kwargs of contourf
    # NOTE: It doens't make sense to use functools.partial here as
    #       contourf is a memberfunction of ax
    # NOTE: zorder sets the rasterization
    #       http://stackoverflow.com/questions/37020842/reducing-size-of-vectorized-contourplot
    cfKwargs = {"cmap"   : cmap  ,\
                "vmax"   : varMax,\
                "vmin"   : varMin,\
                "levels" : levels,\
                "zorder" : -20   ,\
               }

# FIXME: Make these functions: In ax, out ax
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

# FIXME: Only in perpAndPar
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

# FIXME: YOU ARE HERE

# FIXME: Put this on top of decorator
# ??? What are these doing here?
    # Title preparation
    self.helper.rhoTxtDict["value"] =\
            plotNumberFormatter(self.helper.rho[self._xSlice], None)
    self.helper.zTxtDict["value"] =\
            plotNumberFormatter(self.helper.z[self._ySlice], None)
    self.helper.tTxtDict["value"] =\
            plotNumberFormatter(self.helper.t[tInd], None, precision=4)

    # Titles
    polTitle =\
        self.helper.rhoTxtDict["constRhoTxt"].\
            format(self.helper.rhoTxtDict)
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
















#{{{updateColorbar
def updateColorbar(fig, cBarPlane, cBarAx, fluct):
    #{{{docstring
    """
    Updates the colorbar.

    Parameters
    ----------
    fig : Figure
        The figure where the colorbar belongs.
    cBarPlane : Axis
        The axis where the colorbar reads the data from.
    cBarAx : Axis
        The axis where the colorbar belongs.
    fluct: bool
        Whether or not the fluctuations are being plotted.
    """
    #}}}

    # Set the colorbar
    # Clear the axis
    # http://stackoverflow.com/questions/39472017/how-to-animate-the-colorbar-in-matplotlib/39596853
    cBarAx.cla()
    if fluct:
        # Create the ticks (11 with 0 in the center)
        nTicks = 11
        ticks  = np.linspace(self._varMin, self._varMax, nTicks)
        # Enforce the center one to be 0 (without round off)
        ticks[int((nTicks - 1)/2)] = 0
    else:
        ticks = None
    fig.colorbar(cbarPlane,\
                 cax    = self._cBarAx,\
                 ticks  = ticks,\
                 format = FuncFormatter(plotNumberFormatter))
#}}}
