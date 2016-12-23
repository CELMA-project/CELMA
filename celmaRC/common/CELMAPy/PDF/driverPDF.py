#!/usr/bin/env python

"""
Contains drivers for the PDF
"""

from ..commonDrivers import CommonPostProcessingDriver
from ..superClasses import PointsSuperClass
from .plotPDF import PlotPDF
from .calcPDF import calcPDF

#{{{DriverPDF
class DriverPDF(PointsSuperClass, CommonPostProcessingDriver):
    """
    Class which handles the PDF data.
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
        self._PDF = None
    #}}}

    #{{{getPDF
    def getPDF(self):
        """Obtain the PDF"""
        # Create the probes
        self._PDF = calcPDF(self._paths,\
                            self._varName,\
                            self._xInd,\
                            self._yInd,\
                            self._zInd,\
                            converToPhysical = self.convertToPhysical,\
                            mode             = self._mode,\
                            tSlice           = self._tSlice,\
                            )
    #}}}

    #{{{plotPDF
    def plotPDF(self):
        """Plots the PDF"""

        # Calculate the probes if not already done
        if self._PDF == None:
            self.getPDF()

        # Create the energyPlotter
        PDFPlotter = PlotPDF(\
                self._paths                                ,\
                self._PDF                                  ,\
                convertToPhysical = self.convertToPhysical,\
                showPlot          = self._showPlot         ,\
                savePlot          = self._savePlot         ,\
                extension         = self._extension        ,\
                savePath          = self._savePath         ,\
                pltSize           = self._pltSize          ,\
                                  )

        PDFPlotter.plotPDF()
    #}}}
#}}}
