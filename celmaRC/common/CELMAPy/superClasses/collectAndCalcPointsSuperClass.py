#!/usr/bin/env python

"""
Contains parent class for drivers dealing with single points
measurements
"""

from ..collectAndCalcHelpers import (getUniformSpacing,\
                                     getEvenlySpacedIndices,\
                                     findLargestRadialGrad)
from boututils.datafile import DataFile
from boutdata import collect
import numpy as np
import os

#{{{CollectAndCalcPointsSuperClass
class CollectAndCalcPointsSuperClass(object):
    """
    Provides a common constructor interface for driver classes which are
    obtaining temporal information in one point.
    """

    #{{{Constructor
    def __init__(self                     ,\
                 paths                    ,\
                 mode            = "fluct",\
                 nPoints         = None   ,\
                 steadyStatePath = None   ,\
                 ):
        #{{{docstring
        """
        This constructor:

        Set up and stores the xInd, yInd, zInd, tInd and mode

        Parameters
        ----------
        paths  : tuple of strings
            The paths to collect from
        mode : str
            Mode to use when collecting the points,
        nPoints : int
            Size of the sequence. Ignored if xInd, yInd and zInd are all
            sequences.
            See xInd for details.
        steadyStatePath: str
            Path to find the gradient in "n"
        """
        #}}}

        # Set member data
        self._paths           = paths
        self._varName         = varName
        self._mode            = mode
        self._steadyStatePath = steadyStatePath
        self._notCalled       = ["setInd", "setVarName"]
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

    #{{{setInd
    def setInd(self, xInd, yInd, zInd, equallySpace = "x", tSlice=None):
        #{{{docstring
        """
        Sets the slices

        Parameters
        ----------
        xInd : [None|int|sequence of ints]
            xInd to use when collecting.
                * If None: This constructor will use the index of the
                           largest gradient in "n" from the steady state
                           path. This value will be the center-index of
                           nPoints with equidistant spacing
                * If int: If yInd and zInd are given as sequences:
                          The same as None, but the center index will be
                          given by the input int.
                          Else: The value will be copied nPoints times
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
            if type(nPoints) != int:
                message="nPoints will be used, but is specified with a {} type"
                raise ValueError(message.format(type(nPoints)))

        if type(xInd) == int:
            if equallySpace == "x":
                xInd = getEvenlySpacedIndices(paths[0], "x", xInd, nPoints)
            else:
                xInd = (xInd,)*nPoints
        if type(yInd) == int:
            if equallySpace == "y":
                yInd = getEvenlySpacedIndices(paths[0], "y", yInd, nPoints)
            else:
                yInd = (yInd,)*nPoints
        if type(zInd) == int:
            if equallySpace == "z":
                zInd = getEvenlySpacedIndices(paths[0], "z", zInd, nPoints)
            else:
                zInd = (zInd,)*nPoints

        # Guard
        if (len(xInd) != len(yInd)) or\
           (len(xInd) != len(zInd)) or\
           (len(yInd) != len(zInd)):
            raise ValueError("Mismatch in dimension of xInd, yInd and zInd")

        if tSlice is not None:
            if type(tSlice) == int:
                if type(nPoints) != int:
                    message=("nPoints has the wrong type and is needed "
                             "for tSlice multiplication")
                raise ValueError(message)

                tSlice = (tSlice,)*nPoints
            else:
                if len(xInd) != len(tSlice):
                    raise ValueError("Mismatch in dimension of tInd and xInd, yInd and zInd")
    #}}}
#}}}
