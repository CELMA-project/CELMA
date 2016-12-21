#!/usr/bin/env python

"""
Contains class for collecting and calculating the 2D fields
"""

from ..superClasses import CollectAndCalcFieldsSuperClass
from ..collectAndCalcHelpers import (addLastThetaSlice,\
                                     collectiveCollect,\
                                     collectTime,\
                                     get2DMesh,\
                                     polAvg,\
                                     slicesToIndices,\
                                     calcN,\
                                     calcUIPar,\
                                     calcUEPar)
import numpy as np

#{{{CollectAndCalcFields2D
class CollectAndCalcFields2D(CollectAndCalcFieldsSuperClass):
    """
    Class for collecting and calcuating 2D fields
    """

    #{{{constructor
    def __init__(self           ,\
                 *args          ,\
                 fluct  = False ,\
                 mode   = "perp",\
                 **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor
            * Set member data

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        fluct : bool
            If False the raw data is given as an output.
            If True the poloidal average will be subtracted.
        mode : ["perp"|"par"|"pol"|]
            * "perp" - The output field is sliced along a specific z value
            * "par"  - The output field is sliced along a specific theta value
            * "pol"  - The output field is sliced along a specific rho value
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}
        # Guard
        implemented = ("perp", "par", "pol")
        if not(mode in implemented):
            message = "mode '{}' not implemented".format(mode)
            raise NotImplementedError(message)

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._fluct = fluct
        self._mode  = mode
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
                * var    - A 3d array (a 2d spatial array of each time)
                           of the collected variable.
                * varPPi - The field at pi away from the varName field
                           (Only if mode is "par")
                * "X"    - The cartesian x mesh to the field
                * "Y"    - The cartesian Y mesh to the field
                * "time" - The time trace
                * pos    - The position of the fixed index
        """
        #}}}

        # Guard
        if len(self._notCalled) > 0:
            message = "The following functions were not called:\n{}".\
                        format("\n".join(self._notCalled))
            raise RuntimeError(message)

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
            "xInd"         : xInd        ,\
            "yInd"         : yInd        ,\
            "zInd"         : zInd        ,\
                }

        # Collect
        var, time, varPPi =\
            self._collectWrapper(collectKwargs)

        if self.convertToPhysical:
            var  = self.uc.physicalConversion(var , self._varName)
            time = self.uc.physicalConversion(time, "t")
            if self._mode == "par":
                varPPi  = self.uc.physicalConversion(varPPi, self._varName)

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

        If the varName is n, uIPar or uEPar, calculation will be done
        through _calcNonSolvedVars

        Parameters
        ----------
        collectKwargs : dict
           Keyword arguments to use in the collect

        Returns
        -------
        var : array-4d
            The collected array.
        time : array
            The time array
        varPPi : [None|array]
            If mode == "par":
            The collected array pi from var
        """
        #}}}
        if fluct:
            collectKwargs.update({"zInd":None})

        if not(self._varName == "n" or\
               self._varName == "uIPar" or\
               self._varName == "uEPar"):
            var = collectiveCollect(\
                        self._paths      ,\
                        (self._varName, ),\
                        **collectKwargs)
        else:
            var = _calcNonSolvedVars(collectKwargs, not(self.convertToPhysical))

        var  = var[self._varName]
        time = collectTime(self._paths, collectKwargs["tInd"])

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
                var  = var [::self._tSlice.step]
                time = time[::self._tSlice.step]

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

    #{{{_calcNonSolvedVars
    def _calcNonSolvedVars(self, collectKwargs, normalized):
        #{{{docstring
        """
        Calculates variables wich are not solved in the simulation

        Parameters
        ----------
        collectKwargs : dict
            Keyword arguments for collectiveCollect
        normalized : bool
            Whether or not the variables are normalized

        Returns
        -------
        varDict : dict
            Dictionary where the varName variable have been calculated
        """
        #}}}

        lnN = collectiveCollect(self._paths   ,\
                                ("lnN", )     ,\
                                **collectKwargs)
        n = calcN(lnN, normalized)
        if self._varName == "n":
            return n
        else:
            momDensPar = collectiveCollect(self._paths     ,\
                                           ("momDensPar", ),\
                                           **collectKwargs)
            # Calculating uIPar
            uIPar = calcUIPar(momDensPar, n)
            if self._varName == "uIPar":
                return uIPar
            else:
                jPar = collectiveCollect(self._paths   ,\
                                         ("jPar", )    ,\
                                         **collectKwargs)
                uEPar = calcUEPar(uIPar,\
                                  jPar ,\
                                  n    ,\
                                  normalized)
                return uEPar
    #}}}
#}}}
