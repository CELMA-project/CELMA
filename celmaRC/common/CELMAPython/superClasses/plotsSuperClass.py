#!/usr/bin/env python

"""Contains parent class for plotter classes"""

from ..plotHelpers import plotNumberFormatter, seqCMap2, seqCMap3
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.gridspec import GridSpec
import numpy as np
import os

#{{{PlotsSuperClass
class PlotsSuperClass(object):
    """
    Provides a common constructor interface for plotter classes.
    """

    #{{{__init___
    def __init__(self                     ,\
                 path                     ,\
                 convertToPhysical = False,\
                 showPlot          = False,\
                 savePlot          = False,\
                 extension         = "png",\
                 savePath          = "."  ,\
                 pltSize           = None ,\
                 ):
        #{{{docstring
        """
        The constructor for the PlotEnergy object.

        Sets the member data.

        Parameters
        ----------
        path : str
            Paths to create the UnitsConverter object from
        convertToPhysical : bool
            Whether or not to convert to physical units
        showPlot : bool
            If the plots should be displayed.
        savePlot : bool
            If the plots should be saved.
        extension : str
            Extension to use on the plots
        savePath : str
            Path to save destination. Must exist.
        pltSize : tuple
            Size of the plots given as (x, y)
        """
        #}}}

        # Set the member data
        self._showPlot  = showPlot
        self._savePlot  = savePlot
        self._extension = extension
        self._savePath  = savePath
        self._pltSize   = pltSize

        # Make the UnitsConverter object
        self.uc = UnitsConverter(path, convertToPhysical)
        # Reset convertToPhysical from the results of the uc constructor
        self.convertToPhysical = uc.convertToPhysical

        # Make the PlotHelpers object
        self.ph = PlotHelper(self.convertToPhysical)
        self.ph.makeDimensionStringsDicts(self.uc)
    #}}}
#}}}
