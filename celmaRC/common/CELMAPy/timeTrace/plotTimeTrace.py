#!/usr/bin/env python

"""Class for timeTrace plot"""

from ..superClasses import plotSuperClass
from ..plotHelpers import seqCMap3
import numpy as np
import matplotlib.pyplot as plt
import os

#{{{PlotTimeTrace
class PlotTimeTrace(PlotSuperClass):
    """
    Class which contains the time trace data and the plotting configuration.
    """

    #{{{constructor
    def __init__(self, *args, pltSize = (20,15), **kwargs):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets the color
        * Sets the time label

        Parameters
        ----------
        pltSize : tuple
            The size of the plot
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._colors = seqCMap3(np.linspace(0, 1, len(timeTraces.keys())))

        # Set the time label
        self._timeLabel = self._ph.tTxtDict["tTxtLabel"].\
                          format(self.uc.conversionDict["t"])
    #}}}

    #{{{setData
    def setData(self, timeTraces, mode):
        #{{{docstring
        """
        Sets the time traces to be plotted.

        This function also sets the variable labels and the save name.

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
        self._mode = mode

        # Obtain the varname
        ind  = timeTraces.keys()[0]
        keys = timeTraces[ind].keys()
        self._varName = tuple(var for var in keys if var != "time")[0]

        # Set the var label
        pltVarName = self._ph.getVarPltName(self._varname)
        if self._mode == "normal":
            self._varLabel = r"${}$ $[{}]$".\
                    format(pltVarName,\
                           self.uc.conversionDict[self._varName])
            fluctName = ""
        elif self._mode == "fluct":
            self._varLabel = r"$\tilde{{{}}}$ $[{}]$".\
                    format(pltVarName,\
                           self.uc.conversionDict[self._varName])
            fluctName = "fluct"
        else:
            raise NotImplementedError("'{}'-mode not implemented.")

        # Set the fileName
        self._fileName =\
            os.path.join(self._savePath,\
                "{}-{}-{}".format(self._varName, "timeTraces", fluctName))
    #}}}

    #{{{plotSaveShowTimeTrace
    def plotSaveShowTimeTrace(self):
        """
        Performs the actual plotting.
        """

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
        self._ph.makePlotPretty(ax, rotation = 45)

        if self._showPlot:
            plt.show()

        if self._savePlot:
            self._ph.savePlot(fig, self._fileName, (self._leg,))

        plt.close(fig)
    #}}}
#}}}
