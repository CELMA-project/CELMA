#!/usr/bin/env python

"""Class for timeTrace plot"""

from ..superClasses import CollectAndCalcSuperClass
from ..plotHelpers import seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os


#        # Make the UnitsConverter object
#        self.uc = UnitsConverter(path, convertToPhysical)
#        # Reset convertToPhysical from the results of the uc constructor
#        self.convertToPhysical = self.uc.convertToPhysical
#
#        # Make the PlotHelpers object
#        self.ph = PlotHelper(self.convertToPhysical)
#        self.ph.makeDimensionStringsDicts(self.uc)
#{{{PlotTimeTrace
class PlotTimeTrace():
    """
    Class which contains the time trace data and the plotting configuration.
    """

    #{{{__init___
    def __init__(self      ,\
                 timeTraces,\
                 mode      ):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets the member data
        * Prepares the labels
        * Sets the file name

        Parameters
        ----------
        timeTraces : dict
            Dictionary where the keys are on the form "rho,theta,z".
            The value is a dict containing of
            {varName:timeTrace, "time":time}
        mode : ["normal"|"fluct"]
            What mode the input is given in.
        """
        #}}}

        # Set the member data
        self._timeTraces = timeTraces
        self._colors = seqCMap3(np.linspace(0, 1, len(timeTraces.keys())))

        # Obtain the varname
        ind  = timeTraces.keys()[0]
        keys = timeTraces[ind].keys()
        self._varName = tuple(var for var in keys if var != "time")[0]

        # Set the collect path, units converter and plot helper
        self.setCollectPaths(collectPaths, convertToPhysical, uc)

# FIXME: convert?
# FIXME: field1d etc
        # Set the labels
        self._timeLabel = self.ph.tTxtDict["tTxtLabel"].\
                          format(self.uc.conversionDict["t"])

        pltVarName = self.ph.getVarPltName(self._varname)
        if mode == "normal":
            self._varLabel = r"${}$ $[{}]$".\
                    format(pltVarName,\
                           self.uc.conversionDict[self._varName])
        elif mode == "fluct":
            self._varLabel = r"$\tilde{{{}}}$ $[{}]$".\
                    format(pltVarName,\
                           self.uc.conversionDict[self._varName])
        else:
            raise NotImplementedError("'{}'-mode not implemented.")

        self._fileName = "{}.{}".\
            format(os.path.join(self._savePath, "timeTraces"),\
                   self._extension)
    #}}}

    #{{{plotTimeTrace
    def plotTimeTrace(self):
        """ Plots the time traces."""

        # Create the plot
        fig = plt.figure(figsize = self._pltSize)
        ax  = fig.add_subplot(111)

        keys = sorted(self._timeTraces.keys())

        for key, color in keys, self._colors:
            # Make the label
            rho, theta, z = key.split(",")
            label = (r"$\rho={}$ $\theta={}$ $z={}$").format(rho, theta, z)

            ax.plot(self._timeTraces[key]["time"],\
                    self._timeTraces[key][self._varName],\
                    color=color, label=label)

        # Set axis labels
        ax.set_xlabel(self._timeLabel)
        ax.set_ylabel(self._varLabel)

        # Make the plot look nice
        self.ph.makePlotPretty(ax, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self.ph.savePlot(fig, self._fileName, (self._leg,))

        plt.close(fig)
    #}}}
#}}}
