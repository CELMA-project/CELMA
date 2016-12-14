#!/usr/bin/env python

"""
Contains drivers for the probes
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import collectEnergy, PlotEnergy

# Find the max gradient of the variable (subtracts the guard cells)
_, maxGradInd =\
    findLargestRadialGrad(\
      self._varSteadyState[0:1, :, self.yInd:self.yInd+1, 0:1],\
      dx,\
      self._MXG)


# Get radial indices -> equidistanced indices
xind
yind
zind

OR

nXind + center
yind
zind

OR

xind
nYind
zind

# So nYind has higher precedence than yind
# xind etc should still have dimension
# they are called from the outside (e.g. generic driver)
# Make a super class of point type classes, which only contains the init
#--------------




#{{{DriverTimeTrace
class DriverTimeTrace(PointsSuperClass):
    """
    Class which handles the stats and signal data.
    """

    #{{{Constructor
    def __init__(self,\
                 *args,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets pltSize

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._pltSize = (12, 9)

        # Placeholder for the timeTrace
        self._timeTrace = None
    #}}}

    #{{{getTimeTraces
    def getTimeTraces(self):
        """Obtain the timeTrace"""
        # Create the probes
        self._timeTrace = calcTimeTrace(self._paths,\
                                        self._varName,\
                                        self._xInd,\
                                        self._yInd,\
                                        self._zInd,\
                                        self._tSlice,\
                                        self._convertToPhysical)
    #}}}

    #{{{plotTimeTrace
    def plotTimeTrace(self):
        """Plots the timeTrace"""

        # Calculate the probes if not already done
        if self._timeTrace == None:
            self.getTimeTraces()

        # Create the energyPlotter
        energyPlotter = PlotEnergy(\
                self._paths                                ,\
                self._energy                               ,\
                convertToPhysical = self._convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        energyPlotter.plotKinEnergy("ions")
        energyPlotter.plotKinEnergy("electrons")
        energyPlotter.plotPotEnergy()
    #}}}
#}}}
