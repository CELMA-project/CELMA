#!/usr/bin/env python

"""
Contains parent class for drivers dealing with single points
measurements
"""

from .statsAndSignalsDriver import StatsAndSignalsDrivers
from ..statsAndSignals import collectEnergy, PlotEnergy

#{{{PointsSuperClass
class PointsSuperClass(CommonPostProcessingDriver):
    """
    Provides a common constructor interface for driver classes which are
    obtaining temporal information in one point.
    """

    #{{{Constructor
    def __init__(self             ,\
                 *args            ,\
                 paths            ,\
                 tSlice           ,\
                 xInd             ,\
                 yInd             ,\
                 zInd             ,\
                 nPints  = None   ,\
                 varName = "n"    ,\
                 mode    = "fluct",\
                 **kwargs):
        #{{{docstring
        """
        This constructor:

        1. Calls the parent constructor
        2. Sets path and pltSize

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        paths  : tuple of strings
            The paths to collect from
        xInd : [int|sequence of ints]
            xInd to use when collecting.
            If given as int, nPoints must be specified. xInd and nPoints
            will then create a tuple of n equidistant points centred
            around xInd.
            A sequence can also be given directly. In that case nPoints
            will be ignored.
            In both cases, the resulting length of the tuple must match
            the other physical dimensions.
        yInd : [int|sequence of ints]
            The same as xInd, but for the y-index.
        zInd : [int|sequence of ints]
            The same as xInd, but for the z-index.
        nPoints : int
            Number of points in the coordinate given just as an int.
            See xInd for details.
        tSlice : [None|sequence of slices]
            If given this is the  slice of t to use when collecting.
            The length of the sequence must match the other input
            dimensions.
        varName : str
            Name to collect.
        timeTrace : array
            Array of the time trace.
            Will be collected if not given.
        mode : ["normal"|"fluct"]
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._paths   = paths
        self._xInd    = xInd
        self._yInd    = yInd
        self._zInd    = zInd
        self._tSlice  = tSlice
        self._varName = varName
        self._mode    = mode

        # Guard
        if (type(self._xInd) == int and type(self._yInd) == int) or\
           (type(self._xInd) == int and type(self._zInd) == int) or\
           (type(self._yInd) == int and type(self._zInd) == int):
            message = "Only one coordinate can be given as an int at the time"
            raise ValueError(message)

        if type(self._xInd) == int:
            self._xInd = getEvenlySpacedIndices(path, "x", self._xInd, nPoints)
        elif type(self._yInd) == int:
            self._yInd = getEvenlySpacedIndices(path, "y", self._yInd, nPoints)
        elif type(self._zInd) == int:
            self._zInd = getEvenlySpacedIndices(path, "z", self._zInd, nPoints)
    #}}}
#}}}
