#!/usr/bin/env python

"""
Contains classes for plotting the fields
"""

from ..plotHelpers import (PlotHelper,\
                           plotNumberFormatter,\
                           seqCMap,\
                           divCMap,\
                           findLargestRadialGrad)
from ..statsAndSignals import polAvg
from .cylinderMesh import CylinderMesh
from matplotlib import get_backend
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..plotHelpers import safeCollect
import numpy as np
import os
import warnings

#{{{class Plot
class Plot(object):
    """
    Parent class for plotting the results of the CELMA code.

    Handles:

    * Setting the collect options
    * Formation of ticks
    * Preparation of xlabel, ylabel and title
    """

    #{{{Constructor
    def __init__(self                             ,\
                 path                             ,\
                 xguards           = False        ,\
                 yguards           = False        ,\
                 xSlice            = slice(0,None),\
                 ySlice            = slice(0,None),\
                 zSlice            = slice(0,None),\
                 tSlice            = None         ,\
                 convertToPhysical = False        ,\
                 fluctuation         = False        ,\
                 showPlot          = False        ,\
                 savePlot          = True         ,\
                 saveFolder        = None         ,\
                 extension         = "png"        ,\
                 writer            = "ffmpeg"     ,\
                ):
        #{{{docstring
        """
        The constructor:

        1. Sets the plot style
        2. Calcuates rho, theta and z
        3. Collects the time
        4. Collects normalizing parameters if set

        Parameters
        ----------
        path : str
            The path to collect from.
        xguards : bool
            If xguards should be included when collecting.
        yguards : bool
            If yguards should be included when collecting.
        xSlice : slice
            How the data will be sliced in x.
        ySlice : slice
            How the data will be sliced in y.
        zSlice : slice
            How the data will be sliced in z.
        tSlice : slice
            How the data will be sliced in t.
        maxGradRhoFolder : [None | str]
            If this is set, the xSlice will be replaced by the index of
            the largest gradient in rho direction.
        convertToPhysical : bool
            If the physical or normalized units should be plotted.
        fluctuation : bool
            Whether or not the poloidal average should be.
            subtracted from the data.
        showPlot : bool
            If the plot should be displayed.
        savePlot : bool
            If plot should be saved.
        saveFolder : str
            Name of the folder to save plots in.
        extension : str
            Extension of the plot (if the animation is not used).
        writer : str
            Writer to use if the plots are animated.
        """
        #}}}

        # Set the bitrates, fps and codec for ffmpeg (currently magic numbers)
        self._bitrate = -1
        self._fps     = 10
        self._codec   = "h264"

        # Set member data from input
        self._path       = path
        self._xguards    = xguards
        self._yguards    = yguards
        self._showPlot   = showPlot
        self._savePlot   = savePlot
        self._saveFolder = saveFolder
        self._fluctuation  = fluctuation
        self._extension  = extension
        self._writer     = writer

        if writer == "imagemagick":
            self._animExtension = ".gif"
        elif writer == "ffmpeg":
            self._animExtension = ".mp4"
        else:
            message = "Writer {} not implemented".format(writer)
            raise NotImplementedError(message)

        # Get proper indices
        self._xind = self._getIndices(xSlice, "x")
        self._yind = self._getIndices(ySlice, "y")
        self._zind = self._getIndices(zSlice, "z")
        self._tind = self._getIndices(tSlice, "t")

        if type(xSlice) == slice:
            if (xSlice.step is not None):
                message = ("{0}{1}WARNING: xSlice.step not implemented.\n"
                           "Setting to None{1}{0}".format("\n"*3, "!"*3))
                print(message)
                xSlice.step = None
        if type(ySlice) == slice:
            if (ySlice.step is not None):
                message = ("{0}{1}WARNING: ySlice.step not implemented.\n"
                           "Setting to None{1}{0}".format("\n"*3, "!"*3))
                print(message)
                ySlice.step = None
        if type(zSlice) == slice:
            if (zSlice.step is not None):
                message = ("{0}{1}WARNING: ySlice.step not implemented.\n"
                           "Setting to None{1}{0}".format("\n"*3, "!"*3))
                print(message)
                zSlice.step = None

        # Used if we are taking poloidal averages
        self._xSlice = xSlice
        self._ySlice = ySlice
        self._zSlice = zSlice
        # Used to chop the data
        self._tSlice = tSlice

        # Get the time
        t=safeCollect("t_array", path=self._path, tind=self._tind, info=False)

        # Slice in t
        if self._tSlice is not None:
            if self._tSlice.step is not None:
                t = t[::self._tSlice.step]

        # Set frames
        self._frames = len(t)

        # Make the PlotHelper object
        # Public as used in the driver
        self.helper = PlotHelper(path                                 ,\
                                 # Copying as we do not want common memory
                                 t                 = t.copy()         ,\
                                 xguards           = xguards          ,\
                                 yguards           = yguards          ,\
                                 convertToPhysical = convertToPhysical,\
                                 )

        # Set colormap
        if self._fluctuation:
            self._cmap = divCMap
        else:
            self._cmap = seqCMap
    #}}}

    #{{{ _getIndices
    def _getIndices(self, curSlice, dimension):
        """
        Return the slice such that it can be given as an input to "collect"

        Parameters
        ----------
        curSlice : [slice | int | None]
            Current slice to use
        dimension : ["x" | "y" | "z" | "t"]
            The dimension to slice in

        Returns
        -------
        curIndices : sequence (not str)
            A sequence of the start and the stop values in the slice
        """

        if type(curSlice) == slice:
            curIndices = []
            curIndices.append(curSlice.start)
            if curSlice.stop == None:
                # Find the last index
                if dimension == "x" or dimension == "y":
                    dx = safeCollect("dx",\
                                 path=self._path, xguards = self._xguards,\
                                 info=False)
                    dimLen = dx.shape[0]
                if dimension == "y":
                    dy = safeCollect("dy",\
                                 path=self._path, yguards = self._yguards,\
                                 info=False)
                    dimLen = dy.shape[1]
                if dimension == "z":
                    # Subtract 1, as MZ includes the last point
                    dimLen = safeCollect("MZ", path=self._path, info=False) - 1
                if dimension == "t":
                    t = safeCollect("t_array", path=self._path, info=False)
                    dimLen = len(t)
                # Subtract 1 in the end as indices counts from 0
                curIndices.append(dimLen - 1)
            else:
                curIndices.append(curSlice.stop)
        elif curSlice is None:
            curIndices = curSlice
        else:
            curIndices = [curSlice, curSlice]

        # Check for negative indices
        if curIndices is not None:
            for ind in range(len(curIndices)):
                if curIndices[ind] < 0:
                    if dimension == "x" or dimension == "y":
                        dx = safeCollect("dx",\
                                     path=self._path, xguards = self._xguards,\
                                     info=False)
                        dimLen = dx.shape[0]
                    if dimension == "y":
                        dy = safeCollect("dy",\
                                     path=self._path, yguards = self._yguards,\
                                     info=False)
                        dimLen = dy.shape[1]
                    if dimension == "z":
                        # Subtract 1, as MZ includes the last point
                        dimLen =\
                            safeCollect("MZ", path=self._path, info=False) - 1
                    if dimension == "t":
                        t = safeCollect("t_array", path=self._path, info=False)
                        dimLen = len(t)
                    # Subtract 1 in the end as indices counts from 0
                    realInd = dimLen + curIndices[ind] - 1
                    if realInd < 0:
                        message  = ("Index {0} out of range for {1}"
                                    ", as {1} has only {2} elements").\
                            format(curIndices[ind], dimension, dimLen)
                        raise IndexError(message)
                    curIndices[ind] = realInd

        return curIndices
    #}}}
#}}}
