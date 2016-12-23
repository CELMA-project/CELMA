#!/usr/bin/env python

"""
Contains parent class for drivers dealing with single points
measurements
"""

from ..collectAndCalcHelpers import (getUniformSpacing,\
                                     getEvenlySpacedIndices,\
                                     findLargestRadialGradN)
from .driverPostProcessingSuperClass import DriverPostProcessingSuperClass
from boututils.datafile import DataFile
from boutdata import collect
import numpy as np
import os

#{{{CollectAndCalcPointsSuperClass
class CollectAndCalcPointsSuperClass(DriverPostProcessingSuperClass):
    """
    Provides a common constructor interface for driver classes which are
    obtaining temporal information in one point.
    """

    #{{{Constructor
    def __init__(self                     ,\
                 *args                    ,\
                 mode            = "fluct",\
                 nPoints         = None   ,\
                 **kwargs                 ):
        #{{{docstring
        """
        This constructor:

        * Calls the parent constructor
        * Sets the member data

        Parameters
        ----------
        *args : positional arguments
            See the parent constructor for details.
        mode : str
            Mode to use when collecting the points,
        nPoints : int
            Size of the sequence. Ignored if xInd, yInd and zInd are all
            sequences.
            See xInd for details.
        **kwargs : keyword arguments
            See the parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        # Set member data
        self._mode      = mode
        self._nPoints   = nPoints
        self._notCalled = ["setInd", "setVarName"]
    #}}}

    #{{{setVarName
    def setVarName(self, varName):
        #{{{docstring
        """
        Sets the varName

        Parameters
        ----------
        varName : str
            Variable to collect
        """
        #}}}
        self._notCalled.remove("setVarName")
        self._varName = varName
    #}}}

    #{{{setIndices
    def setIndices(self, xInd, yInd, zInd, tSlice=None,\
                   equallySpace = "x", steadyStatePath = None):
        #{{{docstring
        """
        Sets the indices

        Parameters
        ----------
        xInd : [None|int|sequence of ints]
            xInd to use when collecting.
                * If None: This constructor will use the index of the
                           largest gradient in "n" from the steady state
                           path. This value will be the center-index of
                           self._nPoints with equidistant spacing
                * If int: If yInd and zInd are given as sequences:
                          The same as None, but the center index will be
                          given by the input int.
                          Else: The value will be copied self._nPoints times
                * If sequence: All the indices are given
            In all cases, the resulting length of the tuple must match
            the other physical dimensions.
        yInd : [int|sequence of ints]
            The same as xInd (except the None possibility), but for the y-index.
        zInd : [int|sequence of ints]
            The same as xInd (except the None possibility), but for the z-index.
        equallySpace : ["x", "y", "z"]
            If there is any ambiguity of which coordinate to equally
            space around one value, this variable will be used.
            Default is "x".
        tSlice : [None|sequence of slices]
            If given this is the  slice of t to use when collecting.
            The length of the sequence must match the other input
            dimensions.
        steadyStatePath: str
            Path to find the gradient in "n". Only used if xInd is None
        """
        #}}}

        self._notCalled.remove("setInd")

        if xInd is None:
            # Find the x index from the largest gradient in n
            if type(steadyStatePath) != str:
                raise ValueError(("When xInd is None, steadyStatePath "
                                  "must be a string"))

            xInd = findLargestRadialGradN(steadyStatePath)
            if equallySpace != "x":
                print(("{0}Warning: equallySpace was set to {1}, but xInd "
                       "was None. Setting equallySpace to 'x'{0}"
                      ).format("!"*3, equallySpace))

            equallySpace = "x"

        # Guard
        if (type(xInd) == int) or (type(yInd) == int) or (type(zInd) == int):
            if type(self._nPoints) != int:
                message="self._nPoints will be used, but is specified with a {} type"
                raise ValueError(message.format(type(self._nPoints)))

        if type(xInd) == int:
            if equallySpace == "x":
                xInd = getEvenlySpacedIndices(self._collectPaths[0],\
                                              "x", xInd, self._nPoints)
            else:
                xInd = (xInd,)*self._nPoints
        if type(yInd) == int:
            if equallySpace == "y":
                yInd = getEvenlySpacedIndices(self._collectPaths[0],\
                                              "y", yInd, self._nPoints)
            else:
                yInd = (yInd,)*self._nPoints
        if type(zInd) == int:
            if equallySpace == "z":
                zInd = getEvenlySpacedIndices(self._collectPaths[0],\
                                              "z", zInd, self._nPoints)
            else:
                zInd = (zInd,)*self._nPoints

        # Guard
        if (len(xInd) != len(yInd)) or\
           (len(xInd) != len(zInd)) or\
           (len(yInd) != len(zInd)):
            raise ValueError("Mismatch in dimension of xInd, yInd and zInd")

        if tSlice is not None:
            if type(tSlice) == int:
                if type(self._nPoints) != int:
                    message=("self._nPoints has the wrong type and is needed "
                             "for tSlice multiplication")
                raise ValueError(message)

                tSlice = (tSlice,)*self._nPoints
            else:
                if len(xInd) != len(tSlice):
                    raise ValueError("Mismatch in dimension of tInd and xInd, yInd and zInd")

        # Set the member data
        self._xInd   = xInd
        self._yInd   = yInd
        self._zInd   = zInd
        self._tSlice = tSlice
    #}}}
#}}}
