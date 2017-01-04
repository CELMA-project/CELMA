#!/usr/bin/env python

"""
Contains functions for animating the 1D fields
"""

from ..superClasses import PlotAnim1DSuperClass
from ..plotHelpers import plotNumberFormatter
import os

#{{{PlotProfAndGradCompare
class PlotProfAndGradCompare(PlotSuperClass):
    """
    Class which contains the profile and gradient data and the plotting
    configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (15,10), **kwargs):
        #{{{docstring
        """
        Constructor for PlotProfAndGradCompare

        * Calls the parent class
        * Sets the spatial title
        * Sets the xlabel

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

        # Set the member data
        self._pltSize = pltSize
    #}}}

    #{{{setData
    def setData(self, profGrad):
        #{{{docstring
        """
        Sets the profile and gradient to be plotted.

        This function also sets the variable labels, colors and the save name.

        Parameters
        ----------
        profGrad : dict
            Dictionary where the keys are:
                * varName      - Name of the variable under investigation
                * steadyVar    - The profile in the last steady state
                                 time point
                * avgVar       - The averaged profile
                * avgVarStd    - The standard deviation of the average
                * DDXSteadyVar - The radial derivative of the steady
                                 state path
                * DDXAvgVar    - The radial derivative of the average
                * DDXAvgVarStd - The standard deviation of the radial
                                 derivative of the average
                * X            - The radial coordinate
                * zPos         - The fixed z position
        """
        #}}}

        # Set the member data
        self._varName  = profGrad.pop("varName")
        self._rho      = profGrad.pop("X")
        self._zPos     = profGrad.pop("zPos")
        self._profGrad = profGrad

# FIXME: Better colors
        # Obtain the color (pad away brigthest colors)
        pad = 3
        self._colors = seqCMap3(np.linspace(0, 1, 2+pad))

        # Prepare the labels
        self._prepareLabels()






        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varName)





        self._varLabel = self._varLabelTemplate.\
            format(pltVarName, **self.uc.conversionDict[self._varName])

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, "timeTraces", self._fluctName))

        if not(timeAx):
            self._fileName += "Indices"
        if (self._sliced):
            self._fileName += "Sliced"

        if self._extension is None:
            self._extension = "png"

        self._fileName = "{}.{}".format(self._fileName, self._extension)
    #}}}

    #{{{_prepareLabels
    def _prepareLabels(self):
        """
        Prepares the labels for plotting.
        """

        # Set var label template
        if self.uc.convertToPhysical:
            unitsOrNormalization = " $[{units}]$"
            norm                 = ""
            gradUnitsOrNorm      = r" $[{units}\mathrm{{m}}^{{-1}}]$"
        else:

            unitsOrNormalization = r"${normalization}$"
            norm                 = r"${normalization}$"
            # NOTE: The normalization always have a divided by before
            #       the end of the string. Therefore the added rho_s
            #       will appear as a "divided by"
            gradUnitsOrNorm      = r"${normalization}/\rho_s$"

        self._profileLabelTemplate = r"${}$" + unitsOrNormalization
        self._gradLabelTemplate    = r"$\partial_\rho {}$" + gradUnitsOrNorm

# FIXME: Not correct here, see drawing
        self._steadyStateLabel = r"$\mathrm{Steady \quad state}$"
        self._avgLabelTemplate = r"$\lange\langle({0}"+norm+\
                                 r"-\lange\langle {0}"+norm+\
                                 r"\rangle_\theta\rangle_t)^2"+\
                                 r"\rangle_\theta\rangle_t$"

        # Set the spatial part of the title
        self._spatTitle = "{}$,$ {}$,$"
        # Set the x-axis label
# FIXME: Check me
        self._xlabel = self._ph.rhoTxtDict["rhoTxtLabel"]
    #}}}




    #{{{setData
    def setData(self, radialDict, figName):
        #{{{docstring
        """
        Sets the radial data and set up the plotting

        Specifically this function will:
            * Set the data
            * Create the figures and axes (through a function call)
            * Set the colors (through a function call)
            * Set the plot order
            * Update the spatial title

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
        figName : str
            Name of the figure
        plotOrder : [None|sequence of str]
            If given: A sequence of the variable names in the order to
            plot them
        """
        #}}}

        self._X        = radialDict.pop("X")
        self._time     = radialDict.pop("time")
        self._vars     = radialDict

        zPos           = radialDict.pop("zPos")
        thetaPos       = radialDict.pop("thetaPos")

        self._figName  = figName

        # Make axes and colors
        self._createFiguresAndAxes()
        self._setColors()

        # Set the plot order
        if plotOrder:
            self._plotOrder = tuple(plotOrder)
        else:
            self._plotOrder = tuple(sorted(list(self._vars.keys())))

        # Check for "ddt" in the variables
        self._ddtPresent = False

        for var in self._vars:
            if "ddt" in var:
                self._ddtVar = var
                self._ddtPresent = True
                # Remove ddt from the plotOrder
                self._plotOrder = list(self._plotOrder)
                self._plotOrder.remove(var)
                self._plotOrder = tuple(self._plotOrder)
                break

        # Update the title
        self._ph.zTxtDict["value"] = plotNumberFormatter(zPos, None)
        zTxt =\
            self._ph.zTxtDict["constZTxt"].format(self._ph.zTxtDict)
        self._ph.thetaTxtDict["value"] = plotNumberFormatter(thetaPos, None)
        thetaTxt =\
            self._ph.thetaTxtDict["constThetaTxt"].format(self._ph.thetaTxtDict)
        self._spatTitle = self._spatTitle.format(zTxt, thetaTxt)
    #}}}

    #{{{plotAndSaveProfile
    def plotAndSaveProfile(self):
        """
        Performs the actual plotting of the radial plane
        """

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                         "{}-{}-{}".format(self._figName, "radial", "1D"))

        # Initial plot
        self._initialRadialPlot()

        # Call the save and show routine
        self.plotSaveShow(self._fig,\
                          self._fileName,\
                          self._updateRadialAxInTime,\
                          len(self._time))
    #}}}

    #{{{_initialRadialPlot
    def _initialRadialPlot(self):
        #{{{docstring
        """
        Initial radial plot.

        The initial radial plot:
            * Calls the generic initial plot routine
            * Sets the title

        Subsequent animation is done with _updateRadialAxInTime.
        """
        #}}}

        # Initial plot the axes
        self._initialPlot()

        # Set the title
        self._ph.tTxtDict["value"] =\
            plotNumberFormatter(self._time[0], None)
        curTimeTxt = self._ph.tTxtDict["constTTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}{}".format(self._spatTitle, curTimeTxt))
    #}}}

    #{{{_updateRadialAxInTime
    def _updateRadialAxInTime(self, tInd):
        #{{{docstring
        """
        Function which updates the data.

        Parameters
        ----------
        tInd : int
            The current t index.
        """
        #}}}

        for line, key in zip(self._lines, self._plotOrder):
            line.set_data(self._X, self._vars[key][tInd,:])

        # If ddt is present
        if self._ddtPresent:
            # NOTE: ddtLines is one longer than plotOrder
            for line, key in zip(self._ddtLines, self._plotOrder):
                line.set_data(self._X, self._vars[key][tInd,:])

            self._ddtLines[-1].set_data(self._X, self._vars[self._ddtVar][tInd,:])

        # Update the title
        self._ph.tTxtDict["value"] = plotNumberFormatter(self._time[tInd], None)
        curTimeTxt = self._ph.tTxtDict["constTTxt"].format(self._ph.tTxtDict)
        self._fig.suptitle("{}{}".format(self._spatTitle, curTimeTxt))
    #}}}
#}}}

#   # Collect stuff
#   #%%
#   %load_ext autoreload
#   %autoreload 2
#
#   from boutdata import collect
#   import numpy as np
#
#   import os, sys
#   sys.path.append("/home/mmag/CELMA-dev/celma/common/")
#   from CELMAPython.statsAndSignals.averages import timeAvg
#
#   import pickle
#
#   # Must include all zind as we are taking a pol average
#   lnN = collect("lnN", yind=[16,16])
#   t = collect("t_array")
#   n = np.exp(lnN)
#   polSliceN = n[:,16,0,:]
#
#   # nt = [n,t]
#   # with open("/home/mmag/nt", "wb") as f:
#   #     pickle.dump(nt, f)
#
#   from CELMAPython.statsAndSignals.averages import timeAvg, polAvg
#
#   tAvgN, tAvgT = timeAvg(n,t)
#   nAvgTZ = polAvg(tAvgN)
#
#   tLen, xLen, yLen, zLen = n.shape
#
#   nTZFluct = np.zeros(n.shape)
#
#   for t in range(tLen):
#       nTZFluct[t,:,:,:] = n[t,:,:,:] - nAvgTZ[0,:,:,:]
#
#   stdDev = np.sqrt(polAvg(timeAvg(nTZFluct**2.0)))
#
#   avgStd = [nAvgTZ, stdDev]
#   with open("/home/mmag/avgStd.pickle", "wb") as f:
#       pickle.dump(avgStd, f)

#   #%%
#   %load_ext autoreload
#   %autoreload 2
#
#   # Find steady state bg
#   from boutdata import collect
#   import numpy as np
#
#   import pickle
#
#   # Must include all zind as we are taking a pol average
#   lnN = collect("lnN", yind=[16,16])
#   t = collect("t_array")
#   n = np.exp(lnN)
#   steadyStateN = n[0,:,0,0]
#
#   with open("steadyStateN.pickle", "wb") as f:
#       pickle.dump(steadyStateN, f)
#   #%%





# Plot stuff
#%%
import numpy as np
import matplotlib.pylab as plt

import pickle

with open("avgStd.pickle", "rb") as f:
    avg, std = pickle.load(f)
with open("steadyStateN.pickle", "rb") as f:
    steadyStateN = pickle.load(f)

avg = avg[0,:,0,0]
std = std[0,:,0,0]

fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)

# Dummy
rho = range(len(avg))
ax1.plot(rho, avg)
ax1.plot(rho, steadyStateN)

ax1.fill_between(rho,\
                avg+std, avg-std,\
                facecolor="blue", edgecolor="none", alpha=0.5)

# Find the gradient
# FIXME: Need to add the gridspacing there as well, this is given as a
#        positional vararg
gradN       = np.gradient(avg, edge_order=2)
gradNSteady = np.gradient(steadyStateN, edge_order=2)
# From error propagation
gradNStd = np.sqrt((np.abs(gradN)**2)*(std**2))

ax2.plot(rho, gradN)
ax2.plot(rho, gradNSteady)
ax2.fill_between(rho,\
                gradN+gradNStd, gradN-gradNStd,\
                facecolor="blue", edgecolor="none", alpha=0.5)

plt.show()
