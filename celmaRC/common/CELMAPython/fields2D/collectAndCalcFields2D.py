#!/usr/bin/env python

"""
Contains class for collecting and calculating the 2D fields
"""

from ..unitsConverter import UnitsConverter
from ..collectAndCalcHelpers import (DimensionsHelper,\
                                     addLastThetaSlice,\
                                     collectiveCollect,\
                                     get2DMesh,\
                                     slicesToIndices,\
                                     polAvg)
import numpy as np

class CollectAndCalcFields2D(object):
    """
    Class for collecting and calcuating 2D fields
    """

    #{{{constructor
    def __init__(self                      ,\
                 paths                     ,\
                 varName                   ,\
                 xSlice                    ,\
                 ySlice                    ,\
                 zSlice                    ,\
                 tSlice            = None  ,\
                 convertToPhysical = True  ,\
                 fluct             = False ,\
                 uc                = None  ,\
                 dh                = None  ,\
                 mode              = "perp",\
                 xguards           = False ,\
                 yguards           = False ,\
                ):
        #{{{docstring
        """
        Constructor for CollectAndCalcFields2D which sets the member
        data and create the UnitsConverter and DimensionsHelper.

        Parameters
        ----------
        paths : tuple of strings
            The paths to collect from
        varName : str
            Variable to collect
        xSlice : [Slice|int]
            The slice of the rho if the data is to be sliced.
            An integer if "mode" is set to "pol".
        ySlice : [Slice|int]
            The slice of the z if the data is to be sliced.
            An integer if "mode" is set to "perp".
        zSlice : [Slice|int]
            The slice of the z if the data is to be sliced.
            An integer if "mode" is set to "par".
        convertToPhysical : bool
            Whether or not to convert to physical units.
        fluct : bool
            If mode is "normal" the raw data is given as an output.
            If mode is "fluct" the fluctuations are given as an output.
        tSlice : [None|Slice]
            Whether or not to slice the time trace
        uc : [None|UnitsConverter]
            If not given, the function will create the instance itself.
            However, there is a possibility to supply this in order to
            reduce overhead.
        dh : [None|DimensionsHelper]
            If not given, the function will create the instance itself.
            However, there is a possibility to supply this in order to
            reduce overhead.
        mode : ["perp"|"par"|"pol"|]
            * "perp" - The output field is sliced along a specific z value
            * "par"  - The output field is sliced along a specific theta value
            * "pol"  - The output field is sliced along a specific rho value
        xguards : bool
            If the ghost points in x should be collected
        xguards : bool
            If the ghost points in y should be collected
        """
        #}}}

        self._paths             = paths
        self._varName           = varName
        self._xSlice            = xSlice
        self._ySlice            = ySlice
        self._zSlice            = zSlice
        self._tSlice            = tSlice
        self._convertToPhysical = convertToPhysical
        self._fluct             = fluct
        self._mode              = mode
        self._xguards           = xguards
        self._yguards           = yguards

        if uc is None:
            # Create the units convertor object
            uc = UnitsConverter(paths[0], convertToPhysical)
        # Toggle convertToPhysical in case of errors
        convertToPhysical = uc.convertToPhysical

        if dh is None:
            # Create the dimensions helper object
            dh = DimensionsHelper(paths[0], uc)

        self._uc = uc
        self._dh = dh
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates 2D fields

        Returns
        -------
        field2D : dict
            Dictionary with the keys:
                * varName    - A 2D array of the collected variable
                * varNamePPi - The field at pi away from the varName field
                               (Only if mode is "par")
                * "X"        - The cartesian x mesh to the field
                * "Y"        - The cartesian Y mesh to the field
                * "time"     - The time trace
                * pos        - The position of the fixed index
        """
        #}}}

        # Initialize output
        field2D = {}

        # Obtain the mesh
        if self._mode == "perp":
            # X_RT, Y_RT
            X, Y = get2DMesh(rho    = self._dh.rho  ,\
                             theta  = self._dh.theta,\
                             mode   = "RT")
        elif self._mode == "par":
            # X_RZ, Y_RZ
            X, Y = get2DMesh(rho  = self._dh.rho,\
                             z    = self._dh.z  ,\
                             mode = "RZ")
        elif self._mode == "pol":
            # X_ZT, Y_ZT
            X, Y = get2DMesh(theta = self._dh.theta,\
                             z     = self._dh.z,\
                             mode  = "ZT")
        else:
            message = "mode '{}' not implemented".format(self._mode)
            raise NotImplementedError(message)

        # Convert to indices
        xInd = slicesToIndices(self._paths[0], self._xSlice, "x",\
                               xguards=self._xguards)
        yInd = slicesToIndices(self._paths[0], self._ySlice, "y",\
                               yguards=self._yguards)
        zInd = slicesToIndices(self._paths[0], self._zSlice, "z")
        tInd = slicesToIndices(self._paths[0], self._tSlice, "t")

        collectGhost = True if (self._xguards or self._yguards) else False

        # Set keyword arguments
        collectKwargs = {\
            "collectGhost" : collectGhost,\
            "tInd"         : tInd        ,\
            "yInd"         : yInd        ,\
            "xInd"         : xInd        ,\
            "zInd"         : zInd        ,\
                }

        # Collect
        var, time, varPPi =\
            self._collectWrapper(collectKwargs)

        if self._convertToPhysical:
            var  = self._uc.physicalConversion(var , self._varName)
            time = self._uc.physicalConversion(time, "t")
            if self._mode == "par":
                varPPi  = self._uc.physicalConversion(varPPi, self._varName)

        # Store the fields
        field2D["X"   ] = X
        field2D["Y"   ] = Y
        field2D["time"] = time

        if "pol" in self._mode:
            field2D[self._varName ] = var[:, 0, :, :]
            field2D["rhoPos"] = self._dh.rho[xInd[0]]
        if "perp" in self._mode:
            field2D[self._varName] = var[:, :, 0, :]
            field2D["zPos" ] = self._dh.z[yInd[0]]
        if "par" in self._mode:
            field2D[self._varName      ] = var   [:, :, :, 0]
            field2D[self._varName+"PPi"] = varPPi[:, :, :, 0]
            field2D["thetaPos"   ] = self._dh.theta[zInd[0]]

        return field2D
    #}}}

    #{{{_collectWrapper
    def _collectWrapper(self, collectKwargs):
        #{{{docstring
        """
        Collects the variable and the time.

        Parameters
        ----------
        collectKwargs : dict
           Keyword arguments to use in the collect

        Returns
        -------
        var : array
            The collected array
        time : array
            The time array
        varPPi : [None|array]
            If mode == "par":
            The collected array pi from var
        """
        #}}}
        if not(self._fluct):
            varTime = collectiveCollect(\
                            self._paths               ,\
                            (self._varName, "t_array"),\
                            **collectKwargs)

        else:
            collectKwargs.update({"zInd":None})
            varTime = collectiveCollect(\
                                self._paths               ,\
                                (self._varName, "t_array"),\
                                **collectKwargs)

        var  = varTime[self._varName]
        time = varTime["t_array"]

        # Get index
        if self._mode == "par":
            # In this case zSlice is an integer
            zInd = self._zSlice
            # Then theta index corresponding to pi
            piInd = round(var.shape[3]/2)

            if zInd > piInd:
                zPPi = zInd - piInd
            else:
                zPPi = zInd + piInd

        # Collect the negative
        if not(self._fluct) and self._mode == "par":
            collectKwargs.update({"zInd":zPPi})
            varPPi = collectiveCollect(\
                            self._paths      ,\
                            (self._varName, ),\
                            **collectKwargs)

            varPPi = varPPi[self._varName]
        else:
            varPPi = None

        if self._xguards:
            # Remove the inner ghost points from the variable
            var = np.delete(var, (0), axis=1)

        # Slice in t
        if self._tSlice is not None:
            if self._tSlice.step is not None:
                var = var[::self._tSlice.step]

        if self._fluct:
            avg = polAvg(var)
            if self._mode == "par":
                # The negative must have the same average, but not the same
                # fluctuations
                varTmp = (var - avg)[:,:,:,zInd:zInd+1]
                varPPi = (var - avg)[:,:,:,zPPi:zPPi+1]
                var    = varTmp
            else:
                var = (var - avg)

        if self._mode != "par":
            # Add the last theta slice
            var = addLastThetaSlice(var, var.shape[0])

        return var, time, varPPi
    #}}}
