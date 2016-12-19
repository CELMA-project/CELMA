#!/usr/bin/env python

# FIXME: ddt and fix max/min (use max/min routine for this)

# Main

# 2: Create the plot grid
# 3: Do the plot
#    - Can use the max min routine :)...awesome

#{{{  PlotAnim1DRadial
class PlotAnim1DRadial(PlotAnim1DSuperClass):
    """
    Class for 1D radial plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (20,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot1DRadial

        * Calls the parent class
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
        self._radialAx = self._fig.add_subplot(111)
        self._cBarAx = make_axes_locatable(self._radialAx).\
                 append_axes('right', '5%', '5%')
        self._radialAx.grid(True)

        # Set the axis title
        self._axTitle = "{}$,$ {}\n"



# FIXME: After setting the data
        # Set the colors
        colorSpace = np.arange(len(self.lines))
        colors = seqCMap2(np.linspace(0, 1, len(colorSpace)))

# FIXME: A function which sets this
        # Calculate the number of rows
        rows = int(np.ceil(len(self.lines)/self._cols))

        # Create the figure
        fig = plt.figure(figsize = self._pltSize)
        gs  = GridSpec(rows, self._cols)


# FIXME: Move to plotting
        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)


        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self._uc.conversionDict[self._varName])

# FIXME: YOU ARE HERE

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
    #}}}

    #{{{setRadialData
    def setRadialData(self, radialDict, savePath):
        #{{{docstring
        """
        Sets the radial data and label to be plotted

        Parameters
        ----------
        radialDict : dict
            Dictionary of the data containing the following keys
                * var        - The collected variables.
                               A 2d array for each variable
                * "X"        - The abscissa of the variable
                * "time"     - The time trace
                * "zPos"     - The z position
                * "thetaPos" - The theta position
        savePath : str
            Destination to save the plot in.
        """
        #}}}

        self._X        = radialDict.pop("X")
        self._time     = radialDict.pop("time")
        self._zPos     = radialDict.pop("zPos")
        self._thetaPos = radialDict.pop("thetaPos")
        self._vars     = radialDict
        self._savePath = savePath
    #}}}

    #{{{plotAndSaveRadialPlane
    def plotAndSaveRadialPlane(self):
        #{{{docstring
        """
        Performs the actual plotting of the radial plane
        """
        #}}}

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._varName, "radial", "1D"))

        # Initial plot (needed if we would like to save the plot)
        self._updateRadialAxInTime(0)

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updateRadialAxInTime,\
                          len(self._time))
    #}}}

    #{{{_updateRadialAxInTime
    def _updateRadialAxInTime(self, tInd):
        #{{{docstring
        """
        Updates the radial axis.

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
        self._radialAx.cla()

        if self._iterableLevels:
            self._cfKwargs.update({"vmax"   : self._vmax  [tInd],\
                                   "vmin"   : self._vmin  [tInd],\
                                   "levels" : self._levels[tInd],\
                                  })

        # Plot the radial plane
        radialPlane = self._radialAx.\
            contourf(self._X_RT, self._Y_RT, self._Z_RT[tInd, :, :],\
                     **self._cfKwargs)

        # Set rasterization order
        self._radialAx.set_rasterization_zorder(self._axRasterization)
        # Draw the grids
        self._radialAx.grid(b=True)
        # Set x and y labels
        self._radialAx.\
            set_xlabel(self._ph.rhoTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)
        self._radialAx.\
            set_ylabel(self._ph.rhoTxtDict["rhoTxtLabel"],\
                       fontsize = self._labelSize)

        # Update the text
        self._updateRadialPlotTxt(tInd)

        # Update the colorbar
        self._updateColorbar(self._fig, radialPlane, self._cBarAx, tInd)

        # Set title
        self._cBar.set_label(self._varLabel, fontsize = self._labelSize)

        # Set equal axis
        self._radialAx.axis("equal")
    #}}}

    #{{{_updateRadialPlotTxt
    def _updateRadialPlotTxt(self, tInd):
        #{{{docstring
        """
        Updates the radialPlane plot by updating the axis, the title and
        formatting the axes

        Parameters
        ----------
        tInd : int
            The index to plot for
        radialAx : Axis
            The axis to update

        See the docstring of plotRadialPlane for details.
        """
        #}}}

        # Set title
        self._ph.zTxtDict["value"] = plotNumberFormatter(self._constZ, None)
        radialTitle = self._ph.zTxtDict["constZTxt"].format(self._ph.zTxtDict)
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[tInd], None, precision=4)
        timeTitle = self._ph.tTxtDict["tTxt"].format(self._ph.tTxtDict)
        self._radialAx.set_title(self._axTitle.format(radialTitle, timeTitle),\
                               fontsize = self._labelSize)

        # Format axes
        self._ph.makePlotPretty(self._radialAx,\
                                xprune   = "both",\
                                yprune   = "both",\
                                legend   = False,\
                                rotation = 20,\
                                )
    #}}}
#}}}

#{{{Plot1DPar
class Plot1DPar(Plot1DSuperClass):
    """
    Class for 1D parallel plotting.

    Handles plot figure, axes, plot data and decorations.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (20,15), **kwargs):
        #{{{docstring
        """
        Constructor for the Plot1DPar

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
                         "{}-{}-{}".format(self._varName, "par", "1D"))

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
