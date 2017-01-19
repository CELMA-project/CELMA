#!/usr/bin/env python

"""
Contains the total flux calculation
"""

from ..calcVelocities import calcRadialExBConstRho
from ..collectAndCalcHelpers import (DimensionsHelper   ,\
                                     calcUEPar          ,\
                                     calcUIPar          ,\
                                     collectConstZ      ,\
                                     collectConstRho    ,\
                                     collectTime        ,\
                                     getGridSizes       ,\
                                     parallelIntegration,\
                                     poloidalIntegration,\
                                     radialIntegration  ,\
                                     slicesToIndices    ,\
                                    )
from ..unitsConverter import UnitsConverter
import numpy as np

#{{{CollectAndCalcTotalFlux
class CollectAndCalcTotalFlux(object):
    """
    Class for collecting and calcuating the total flux
    """

    #{{{constructor
    def __init__(self                        ,\
                 collectPaths                ,\
                 xInd              = None    ,\
                 yInd              = None    ,\
                 tSlice            = None    ,\
                 mode              = "normal",\
                 convertToPhysical = True    ,\
                 ):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor

        Parameters
        ----------
        collectPaths : tuple
            Tuple from where to collect
        xInd : [None|int]
            How to slice in the radial direction.
            If None, the last inner point is selected.
        yInd : [None|int]
            How to slice in the parallel direction.
            If None, the last inner point is selected.
        yInd : [None|int]
        tSlice : [None|slice]
            How to slice in time.
        mode : ["normal"|"fluct"]
            Whether to look at fluctuations or normal data
        convertToPhysical : bool
            Whether or not to convert to physical
        """
        #}}}

        # Set the member data
        self._collectPaths = collectPaths
        self._mode         = mode
        self._tSlice       = tSlice

        # Set the indices
        if xInd is None:
            # NOTE: Indices counting from 0
            self._xInd = getGridSizes(collectPaths[0], "x") - 1
        else:
            self._xInd  = xInd
        if yInd is None:
            # NOTE: Indices counting from 0
            self._yInd = getGridSizes(collectPaths[0], "y") - 1
        else:
            self._yInd  = yInd
        # Get the units converter
        self.uc = UnitsConverter(self._collectPaths[0], convertToPhysical)
        self.convertToPhysical = self.uc.convertToPhysical
        # Get the dimensions helper
        self._dh = DimensionsHelper(self._collectPaths[0], self.uc)

        # Get the tInd trace
        self._tInd = slicesToIndices(self._collectPaths[0], self._tSlice, "t")
    #}}}

    #{{{executeCollectAndCalc
    def executeCollectAndCalc(self):
        #{{{docstring
        """
        Function which collects and calculates the integrated fluxes.

        This fucntion will:
            1. Collect the velocities
            2. Calculate n from collected lnN
            3. Multiply and integrate
            4. Recast to 1D
            5. Out parallel int flux, perpendicular int flux

        Returns
        -------
        intFluxes : dict
            Dictionary on the form:
            {"parElIntFlux"  : parallelElectronIntegratedFlux,
             "parIonIntFlux" : parallelIonIntegratedFlux,
             "perpIntFlux"   : perpendicularIntegratedFlux,
             "time"          : time,
            }
        """
        #}}}

        # Initialize output
        intFluxes = {}

        # Collect densities
        radialN = np.exp(self._collectAndCalcConstRho("lnN"))
        parN    = np.exp(self._collectAndCalcConstZ("lnN"))

        if self.convertToPhysical:
            radialN = self.uc.physicalConversion(radialN, "n")
            parN    = self.uc.physicalConversion(parN   , "n")

        # Collect the parallel velocities
        parMomDensPar = self._collectAndCalcConstZ("momDensPar")
        parJPar       = self._collectAndCalcConstZ("jPar")
        import pdb; pdb.set_trace()
        parIonVel = calcUIPar(parMomDensPar, parN)
        parElVel  = calcUEPar(parIonVel, parJPar, parN,\
                              not(self._collectAndCalcConstZ))

        # Collect the perpendicular velocities
        radialExB = calcRadialExBConstRho(\
                          self._collectPaths                        ,\
                          self._xInd                                ,\
                          self._tSlice                              ,\
                          mode              = self._mode            ,\
                          convertToPhysical = self.convertToPhysical,\
                          )

        # Collect time
        time = collectTime(self._collectPaths, tInd = self._tInd)
        if self.convertToPhysical:
            time = self.uc.physicalConversion(time ,"t")

        if self._tSlice is not None:
            if type(self._tSlice) == slice:
                if self._tSlice.step is not None:
                    time  = time[::self._tSlice.step]

        # Multiply
        radFluxDens    = radialN*radialExB
        parElFluxDens  = parN*parElVel
        parIonFluxDens = parN*parIonVel

        # Integration multipliers
        rho = self._dh.rho[self._xInd]
        dx = self._dh.dx
        dy = self._dh.dy

        # First integration
        intRadFluxDens    = poloidalIntegration(radFluxDens   , rho)
        intParElFluxDens  = poloidalIntegration(parElFluxDens , rho)
        intParIonFluxDens = poloidalIntegration(parIonFluxDens, rho)

        # Storing
        intFluxes["parElIntFlux"]  = radialIntegration  (intRadFluxDens   , dx)
        intFluxes["parIonIntFlux"] = radialIntegration  (intParElFluxDens , dx)
        intFluxes["perpIntFlux"]   = parallelIntegration(intParIonFluxDens, dy)
        intFluxes["time"]          = time

        return intFluxes
    #}}}

    #{{{_collectAndCalcConstZ
    def _collectAndCalcConstZ(self, varName):
        #{{{docstring
        """
        Collects and transforms a variable for a constant z

        Parameters
        ----------
        varName : str
            The variable to collect.

        Returns
        -------
        var : array-4d
            The collected variable
        time : array-1d
            The corresponding time
        """
        #}}}

        tInd = slicesToIndices(self._collectPaths[0], self._tSlice, "t")

        # Collect phi
        var = collectConstZ(self._collectPaths,\
                            varName           ,\
                            self._yInd        ,\
                            tInd=self._tInd)

        # Slice
        if self._tSlice is not None:
            if type(self._tSlice) == slice:
                if self._tSlice.step is not None:
                    var  = var [::self._tSlice.step]

        # Convert to physical units
        if self.convertToPhysical:
            var = self.uc.physicalConversion(var, varName)

        if self._mode == "fluct":
            var = (var - polAvg(var))

        return var
    #}}}

    #{{{_collectAndCalcConstRho
    def _collectAndCalcConstRho(self, varName):
        #{{{docstring
        """
        Collects and transforms a variable for a constant rho

        Parameters
        ----------
        varName : str
            The variable to collect.

        Returns
        -------
        var : array-4d
            The collected variable
        time : array-1d
            The corresponding time
        """
        #}}}

        # Collect phi
        var = collectConstRho(self._collectPaths,\
                              varName           ,\
                              self._xInd        ,\
                              tInd = self._tInd)

        # Slice
        if self._tSlice is not None:
            if type(self._tSlice) == slice:
                if self._tSlice.step is not None:
                    var  = var [::self._tSlice.step]

        # Convert to physical units
        if self.convertToPhysical:
            var = self.uc.physicalConversion(var, varName)

        if self._mode == "fluct":
            var = (var - polAvg(var))

        return var
    #}}}
#}}}
