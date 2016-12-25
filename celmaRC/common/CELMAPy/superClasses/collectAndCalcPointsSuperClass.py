#!/usr/bin/env python

"""
Contains class for collecting and calculating points
"""

from ..collectAndCalcHelpers import (findLargestRadialGradN,\
                                     getEvenlySpacedIndices)
from .collectAndCalcSuperClass import CollectAndCalcSuperClass

#{{{CollectAndCalcPointsSuperClass
class CollectAndCalcPointsSuperClass(CollectAndCalcSuperClass):
    """
    Super class for field collection.

    Handles setting of data for field collection.
    """

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        Super constructor for collection and calculation of fields.

        This constructor will:
            * Call the parent class
            * Set the member data

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

        self._notCalled.append("setIndices")
    #}}}

    #{{{setIndices
    def setIndices(self, xInd, yInd, zInd, tSlice=None,\
                   nPoints = None, equallySpace = "x", steadyStatePath = None):
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
                           nPoints with equidistant spacing
                * If int: If yInd and zInd are given as sequences:
                          The same as None, but the center index will be
                          given by the input int.
                          Else: The value will be copied nPoints times
                * If sequence: All the indices are given
            In all cases, the resulting length of the tuple must match
            the other physical dimensions.
        yInd : [int|sequence of ints]
            The same as xInd (except the None possibility) for the y-index.
        zInd : [int|sequence of ints]
            The same as xInd (except the None possibility) for the z-index.
        tSlice : [None|sequence of slices]
            If given this is the  slice of t to use when collecting.
            The length of the sequence must match the other input
            dimensions.
        nPoints : [None|int]
            Size of the sequence. Ignored if xInd, yInd and zInd are all
            sequences.
            See xInd for details.
        equallySpace : ["x", "y", "z"]
            If there is any ambiguity of which coordinate to equally
            space around one value, this variable will be used.
            Default is "x".
        steadyStatePath: str
            Path to find the gradient in "n". Only used if xInd is None
        """
        #}}}

        self._notCalled.remove("setIndices")

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
                xInd = getEvenlySpacedIndices(self._collectPaths[0],\
                                              "x", xInd, nPoints)
            else:
                xInd = (xInd,)*nPoints
        if type(yInd) == int:
            if equallySpace == "y":
                yInd = getEvenlySpacedIndices(self._collectPaths[0],\
                                              "y", yInd, nPoints)
            else:
                yInd = (yInd,)*nPoints
        if type(zInd) == int:
            if equallySpace == "z":
                zInd = getEvenlySpacedIndices(self._collectPaths[0],\
                                              "z", zInd, nPoints)
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
                    message = ("Mismatch in dimension of tInd and " 
                               "xInd, yInd and zInd")
                    raise ValueError(message)

        # Set the member data
        self._xInd   = xInd
        self._yInd   = yInd
        self._zInd   = zInd
        self._tSlice = tSlice
    #}}}
#}}}
