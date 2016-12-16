#!/usr/bin/env python

"""
Contains functions for plotting the 2D fields
"""

#{{{plotPerpPlane
def plotPerpPlane(X_RT, Y_RT, Z_RT, time, constZ, ph,\
                  savePath = ".test", \
                  cfKwargs={}, txtKwargs={}):
    #{{{docstring
    """
    Performs the actual plotting of the perpendicular plane

    Parameters
    ----------
    X_RT : array
        A 2d mesh of the Cartesian x coordinates
    Y_RT : array
        A 2d mesh of the Cartesian y coordinates
    Z_RT : array
        A 3d array of yhe vaules for each point in x and y for each time
    ph : PlotHelper
        The plot helper
    savePath : str
        Path (excluding extension) to use if the plot or animation is
        saved.
    show : bool
        Whether or not the plot is to be displayed.
    save : bool
        Whether or not to save the plot.
    cfKwargs : dict
        Extra keywords for the contourf plot
    txtKwargs : dict
        Extra keywords for setting the text on the title
    """
    #}}}

    # Create figure and axes
    pltSize = (20,15)
    fig = plt.figure(figsize = pltSize)
    perpAx = fig.add_subplot(111)
    cBarAx = make_axes_locatable(perpAx).\
             append_axes('right', '5%', '5%')
    perpAx.grid(True)

    # Initial plot (needed if we would like to save the plot)
    updatePerpPlane(0, perpAx, X_RT, Y_RT, Z_RT,\
                    time, constZ, ph, cfKwargs, txtKwargs)

    # Call the save and show routine
    fargs = (perpAx, X_RT, Y_RT, Z_RT, time, constZ, ph, cfKwargs, txtKwargs)
    plotSaveAndShow(fig, updatePerpPlane, fargs, Z_RT.shape[0],\
                    show = show, save = save, savepath = savepath)
#}}}

#{{{updatePerpPlanePlot
def updatePerpPlanePlot(tInd, perpAx, X_RT, Y_RT, Z_RT,\
                        time, constZ, ph, cfKwargs, txtKwargs):
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

    perpPlane = updatePerpAx(perpAx, X_RT, Y_RT, Z_RT[tInd, :, :], ph, cfKwargs)

    # Set title
    ph.zTxtDict["value"] = plotNumberFormatter(constVal, None)
    perpTitle = ph.zTxtDict["constZTxt"].format(ph.zTxtDict)
    ph.tTxtDict["value"] = plotNumberFormatter(t[tInd], None, precision=4)
    timeTitle = ph.tTxtDict["tTxt"].format(ph.tTxtDict)
    perpTxt = perpAx.text(1.5, 1.05,\
                          "{}$,$ {}".format(perpTitle, timeTitle),\
                          transform = self._perpAx.transAxes,\
                          **txtKwargs)

    # Format axes
    perpAx.get_xaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
    perpAx.get_yaxis().set_major_formatter(FuncFormatter(plotNumberFormatter))
    perpAx.axis("equal")

    updateColorbar(fig, perpPlane, cBarAx, fluct)
#}}}

#{{{updatePerpAx
def updatePerpAx(perpAx, X_RT, Y_RT, Z_RT, ph, cfKwargs, labelSize=35):
    #{{{docstring
    """
    Updates the perpoidal axis by cleaning it, and reset the plot

    This also sets the rasterization, the grid and the labels.

    Parameters
    ----------
    perpAx : Axis
        The axis to update
    labelSize : int
        Size of the labelfont

    See the docstring of updatePerpPlanePlot for details.

    Returns
    -------
    perpPlane : contour-like
        The object which can be used to find a colorbar
    """
    #}}}

    # Clear previous axis
    perpAx.cla()

    # Plot the perpoidal plane
    perpPlane = ax.\
        contourf(self._cyl.X_RT, self._cyl.Y_RT, Z_RT, **cfKwargs)

    # Set rasterization order
    perpAx.set_rasterization_zorder(-10)
    # Draw the grids
    perpAx.grid(b=True)
    # Set x and y labels
    perpAx.set_xlabel(ph.zTxtDict["rhoTxtLabel"], fontsize = labelSize)
    perpAx.set_ylabel(ph.zTxtDict["rhoTxtLabel"], fontsize = labelSize)

    return perpPlane
#}}}

#{{{updateParAx
def updateParAx(parAx, X_RZ, Y_RZ, Z_RZ, Z_RZ_PPi,\
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
    contourf(-X_RZ, Y_RZ, Z_RZ_PPi, **cfKwargs)

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
def updatePolAx(polAx, X_ZT, Y_ZT, Z_ZT, ph, cfKwargs={}, labelSize=35):
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
        contourf(self._cyl.X_ZT, self._cyl.Y_ZT, Z_ZT.transpose(), **cfKwargs)

    # Set rasterization order
    polAx.set_rasterization_zorder(-10)
    # Draw the grids
    polAx.grid(b=True)
    # Set x and y labels
    polAx.set_xlabel(r"$\theta$", fontsize = labelSize)
    polAx.set_ylabel(ph.zTxtDict["zTxtLabel"], fontsize = labelSize)

    return polPlane
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
    cBarPlane : contour-like
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
                 cax    = cBarAx,\
                 ticks  = ticks,\
                 format = FuncFormatter(plotNumberFormatter))
#}}}

#{{{
def plotSaveAndShow(fig, func, fargs, frames,\
                    bitrate = -1, fps = 10, codec = "h264",\
                    show = False, save = True, savepath=".test",\
                    extension = None):
    #{{{docstring
    """
    Saves and closes the plot.

    Parameters
    ----------
    fig : Figure
        Figure to save.
    func : function
        The function to use for generating the animation.
    fargs : tuple
        Tuple of the arguments to use in the animation.
    frames : int
        Number of frames.
        If this is less than one, a normal plot will be made.
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
    show : bool
        Whether or not the plot is to be displayed.
    save : bool
        Whether or not to save the plot.
    savePath : str
        Path (excluding extension) to use if the plot or animation is
        saved.
    extension : [None|str]
        Overrides default extension if set.
    """
    #}}}

    if frames > 1:

        # Animate
        anim = animation.FuncAnimation(fig            ,\
                                       func           ,\
                                       fargs  = fargs ,\
                                       frames = frames,\
                                       blit   = False ,\
                                       )

        if savePlot:
            FFMpegWriter = animation.writers['ffmpeg']
            writer = FFMpegWriter(bitrate = bitrate,\
                                  fps     = fps    ,\
                                  codec   = codec)

            if extension is None:
                extension = ".mp4"

            # Save the animation
            anim.save(savepath + extension, writer = writer)
            print("Saved to {}{}".format(savepath, extension))
    else:
        if savePlot:

            if extension is None:
                extension = ".pdf"

            # Save the figure
            fig.savefig("{}.{}".format(savepath, extension),\
                        transparent = True  ,\
                        bbox_inches = "tight",\
                        pad_inches  = 0      ,\
                        )
            print("Saved to {}.{}".format(savepath, extension))

    if showPlot:
        fig.show()

    plt.close(fig)
#}}}
