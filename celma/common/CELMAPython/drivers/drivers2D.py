#!/usr/bin/env python

"""
Contains drivers for plotting 2D plots
"""

from ..modelSpecific import varAndPlotNames
from ..fieldPlotters import Plot2D
from .fieldPlottersDriver import FieldPlottersDriver
import numpy as np
from multiprocessing import Process
from boutdata import collect

#{{{Drivers2D
class Drivers2D(FieldPlottersDriver):
    """
    Class which handles the 2D plots of the fields.
    """

    #{{{Constructor
    def __init__(self                    ,\
                 *args                   ,\
                 varName           = None,\
                 var               = None,\
                 pltName           = None,\
                 varMax            = None,\
                 varMin            = None,\
                 varyMaxMin        = None,\
                 axisEqualParallel = True,\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets additional member data.

        Parameters
        ----------
        *args : positional arguments
            See the constructor of FieldPlottersDriver for details.
        varName : str
            Name of the field which is going to be collected (if var is
            not given).
        var : [None | array]
            The variable to plot.
        pltName : str
            Name of the plot written in LaTeX format, but without the $.
        varMax : float
            Setting a hard upper limit z-axis in the plot.
        varMin : float
            Setting a hard lower limit z-axis in the plot.
        varyMaxMin : bool
            Whether or not the limits of the z-axis should be set
            to the max/min of the current timestep or not.
        **kwargs : keyword arguments
            See the constructor of FieldPlottersDriver for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._varName           = varName
        self._var               = var
        self._pltName           = pltName
        self._varMax            = varMax
        self._varMin            = varMin
        self._varyMaxMin        = varyMaxMin
        self._axisEqualParallel = axisEqualParallel
    #}}}

    #{{{allMainFields2DDriver
    def allMainFields2DDriver(self):
        """
        Driver for all 2D plots.
        """

        if self._useSubProcess:
            #{{{ Function call through subprocess
            for varName, plotName in varAndPlotNames:
                self._varName = varName
                self._pltName = plotName
                Process(\
                        target = self.single2DDriver,\
                        args   = ()  ,\
                        kwargs = {}
                       ).start()
            #}}}
        else:
            #{{{ Normal function call
            # Do the plotting
            for varName, plotName in varAndPlotNames:
                self._varName = varName
                self._pltName = plotName
                self.single2DDriver()
            #}}}
    #}}}

    #{{{single2DDriver
    def single2DDriver(self):
        """
        Driver for a single 2D plot.
        """

        try:
            # Make the plotter object
            plotter = Plot2D(self._path                                 ,\
                             self._varName                              ,\
                             var               = self._var              ,\
                             xguards           = self._xguards          ,\
                             yguards           = self._yguards          ,\
                             xSlice            = self._xSlice           ,\
                             ySlice            = self._ySlice           ,\
                             zSlice            = self._zSlice           ,\
                             tSlice            = self._tSlice           ,\
                             subPolAvg         = self._subPolAvg        ,\
                             convertToPhysical = self._convertToPhysical,\
                             showPlot          = self._showPlot         ,\
                             savePlot          = self._savePlot         ,\
                             saveFolder        = self._saveFolder       ,\
                             varMax            = self._varMax           ,\
                             varMin            = self._varMin           ,\
                             varyMaxMin        = self._varyMaxMin       ,\
                             axisEqualParallel = self._axisEqualParallel,\
                             extension         = self._extension        ,\
                            )
        except (KeyError, ValueError) as collectError:

            # Get the tind
            if self._tSlice is not None:
                self._tind = [self._tSlice.start]
                if self._tSlice.stop == None:
                    t = collect("t_array", path=self._path, info=False)
                    dimLen = len(t)
                    # Subtract 1 in the end as indices counts from 0
                    self._tind.append(dimLen - 1)
            else:
                self._tind = None

            if self._varName == "n":
                #{{{n
                lnN = collect("lnN"                  ,\
                              path    = self._path   ,\
                              yguards = self._yguards,\
                              xguards = self._xguards,\
                              tind    = self._tind   ,\
                              info    = False        ,\
                              )

                self._var = np.exp(lnN)
                #}}}
            if self._varName == "jPar":
                #{{{jPar
                lnN = collect("lnN"                  ,\
                              path    = self._path   ,\
                              yguards = self._yguards,\
                              xguards = self._xguards,\
                              tind    = self._tind   ,\
                              info    = False        ,\
                              )

                uEPar = collect("uEPar"                ,\
                                path    = self._path   ,\
                                yguards = self._yguards,\
                                xguards = self._xguards,\
                                tind    = self._tind   ,\
                                info    = False        ,\
                                )

                uIPar = collect("uIPar"                ,\
                                path    = self._path   ,\
                                yguards = self._yguards,\
                                xguards = self._xguards,\
                                tind    = self._tind   ,\
                                info    = False        ,\
                                )

                self._var = np.exp(lnN)*(uIPar - uEPar)
                #}}}

            # Make the plotter object
            plotter = Plot2D(self._path                                 ,\
                             self._varName                              ,\
                             var               = self._var              ,\
                             xguards           = self._xguards          ,\
                             yguards           = self._yguards          ,\
                             xSlice            = self._xSlice           ,\
                             ySlice            = self._ySlice           ,\
                             zSlice            = self._zSlice           ,\
                             tSlice            = self._tSlice           ,\
                             subPolAvg         = self._subPolAvg        ,\
                             convertToPhysical = self._convertToPhysical,\
                             showPlot          = self._showPlot         ,\
                             savePlot          = self._savePlot         ,\
                             saveFolder        = self._saveFolder       ,\
                             varMax            = self._varMax           ,\
                             varMin            = self._varMin           ,\
                             varyMaxMin        = self._varyMaxMin       ,\
                             axisEqualParallel = self._axisEqualParallel,\
                            )

        plotter.plotDriver(self._pltName, savePath = self._savePath)
    #}}}
#}}}
