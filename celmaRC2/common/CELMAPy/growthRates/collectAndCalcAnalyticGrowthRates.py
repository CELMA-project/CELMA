#!/usr/bin/env python

from ..collectAndCalcHelpers import (DimensionsHelper      ,\
                                     collectSteadyN        ,\
                                     DDX                   ,\
                                     getUniformSpacing     ,\
                                     getScanValue          ,\
                                     findLargestRadialGradN,\
                                    )
from .analyticalGrowthRates import (calcOmCE         ,\
                                    calcOmStar       ,\
                                    calcPecseliB     ,\
                                    calcSigmaPar     ,\
                                    calcUDE          ,\
                                    pecseliAnalytical,\
                                   )
from .unitsConverter import UnitsConverter
from boututils.datafile import DataFile
import pandas as pd
import scipy.constants as cst
import os

#{{{CollectAndCalcAnalyticGrowthRates
class CollectAndCalcAnalyticGrowthRates(object):
    """
    Class for collecting and calculating variables needed for
    calculations of the analytical growth rates.
    """

    #{{{constructor
    def __init__(self, steadyStatePaths, scanParameter, yInd):

        #{{{docstring
        """
        This constructor will:
            * Set the member data

        Parameters
        ----------
        steadyStatePaths : tuple
            Path to the steady state simulation.
            The tuple must be ordered according to the scan order in
            scanCollectPaths.
        scanParameter : str
            String segment representing the scan
        yInd : int
            y-index to collect from
        """
        #}}}

        self._steadyStatePaths = steadyStatePaths
        self._scanParameter    = scanParameter
        self._yInd             = yInd
        self._maxRhoInd        = None
    #}}}

    #{{{_setUcAndDh
    def _setUcAndDh(self, path):
        #{{{docstring
        """
        Sets the unitsConverter and the dimensionsHelper

        Parameters
        ----------
        path : str
            The path to use for the objects
        """
        #}}}
        # Create the units convertor object
        self.uc  = UnitsConverter  (path, True)
        self._dh = DimensionsHelper(path, self.uc)
    #}}}

    #{{{_collectForPecseliSemiAnalytical
    def _collectForPecseliSemiAnalytical(self, path, yInd):
        #{{{docstring
        """
        Collects variables needed for analytical calcualtion of the
        Pecseli dispersion relation.

        See also
        --------
        _calcPecseliSemiAnalytical

        Parameters
        ----------
        path : str
            Path to collect from
        yInd : int
            The parallel index to slice the data

        Returns
        -------
        omCE : float
            The electron cyclotron resonance.
        omCI : float
            The ion cyclotron resonance.
        rhoS : float
            The hybrid radius.
        rhoMax : float
            The position of the maximum gradient.
        n : float
            The density at the maximum gradient.
        dndx : float
            The density gradient at the position of the maximum
            gradient.
        uDE : float
            The electron diamagnetic drift.
        nuEI : float
            The electron-ion collision frequency.
        """
        #}}}
        # Get the position of the largest gradient
        indMaxGrad = findLargestRadialGradN(path)
        rhoMax = self._dh.rho[indMaxGrad]

        # Calculate B
        mi   = self.uc.getNormalizationParameter("mi")
        omCI = self.uc.getNormalizationParameter("omCI")*mi/cst.e
        B    = omCI*mi/cst.e

        # Calculation of omStar
        # NOTE: For now: Duplication of work as gradient in already
        #       calculated in findLargestRadialGradN
        n    = collectSteadyN(path, yInd)
        dx   = getUniformSpacing(path, "x")
        dndx = DDX(n, dx)
        n    = n   [0,indMaxGrad,yInd,0]
        dndx = dndx[0,indMaxGrad,yInd,0]
        Te   = self.uc.getParamFromNormDict("Te0")
        uDE  = calcUDE(Te, B, n, dndx)

        # Needed for calculation of sigmaPar
        omCE = calcOmCE(B)
        with DataFile(os.path.join(path,"BOUT.dmp.0.nc")) as f:
            nuEI = f.read("nuEI")

        # Needed for calculation of b
        rhoS = self.uc.getNormalizationParameter("rhoS")

        return omCE, omCI, rhoS, rhoMax, n, dndx, uDE, nuEI
    #}}}

    #{{{getData
    def getData(self, nModes = 7):
        #{{{docstring
        """
        Makes a DataFrame of the growth rates and angular velocities.

        NOTE:
            * Assumes singly ionization
            * Assumes that kz is twice the simulation box

        Parameters
        ----------
        nModes : int
            Number of modes.

        Returns
        -------
        analyticalGRDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "analyticalGR"
                * "angularVelocity"
            over the observation "modeNr" over the observation "Scan"
        paramDataFrame : DataFrame
            DataFrame consisting of the variables (measured properties):
                * "omCE"
                * "omCI"
                * "rhoS"
                * "rhoMax"
                * "n"
                * "dndx"
                * "uDE"
                * "nuEI"
            over the observation "Scan".
        positionTuple : tuple
            The tuple containing (rho, z).
            Needed in the plotting routine.
        uc : Units Converter
            The units converter used when obtaining the fourier modes.
            Needed in the plotting routine.
        """
        #}}}

        # For ideas on how to append a DataFrame, see:
        # https://github.com/pandas-dev/pandas/blob/master/doc/cheatsheet/Pandas_Cheat_Sheet.pdf
        # http://stackoverflow.com/questions/10715965/add-one-row-in-a-pandas-dataframe

        singleIndexTuple = []
        multiIndexTuple  = []
        fullDict = {"analyticalGR":[], "angularVelocity":[]}

        keys = ("omCE", "omCI", "rhoS", "rhoMax", "n", "dndx", "uDE", "nuEI")
        paramDict = {key:[] for key in keys}

        # Loop over the folders
        for steadyStatePath in self._steadyStatePaths:
            # Set the units converter and the dimension helper
            self._setUcAndDh(steadyStatePath)
            # Obtain the scan value
            scanValue = getScanValue(steadyStatePath, self._scanParameter)

            # Assume kz is double of the cylinder heigth (observed in
            # the fluctuations)
            kz = 2*(self._dh.z[-1] + 0.5*(self._dh.z[-1] - self._dh.z[-2]))

            # Collect variables
            omCE, omCI, rhoS, rhoMax, n, dndx, uDE, nuEI =\
                self._collectForPecseliSemiAnalytical(steadyStatePath,\
                                                      self._yInd)

            # Update the single index tuple, insert into parameter dictionary
            singleIndexTuple.append(scanValue)
            paramDict["omCE"]  .append(omCE)
            paramDict["omCI"]  .append(omCI)
            paramDict["rhoS"]  .append(rhoS)
            paramDict["rhoMax"].append(rhoMax)
            paramDict["n"]     .append(n)
            paramDict["dndx"]  .append(dndx)
            paramDict["uDE"]   .append(uDE)
            paramDict["nuEI"]  .append(nuEI)

            for m in range(len(1, nModes+1)):
                # Fill the multiIndexTuple and the dict
                multiIndexTuple.append((scanValue, m))

                # Calculate the dispersion relation
                ky = m/rhoMax
                omStar = calcOmStar(ky, uDE)
                b = calcPecseliB(ky, rhoS)
                sigmaPar = calcSigmaPar(ky, kz, omCE, omCI, nuEI)
                om = pecseliAnalytical(omStar, b, sigmaPar)

                fullDict["analyticalGR"].append(om.imag)
                fullDict["angularVelocity"].append(om.real)

        # Make the data frames
        paramDataFrame = pd.DataFrame(paramDict, index=singleIndexTuple)

        analyticalGRDataFrame =\
            pd.DataFrame(fullDict,\
                         index=pd.MultiIndex.from_tuples(\
                            multiIndexTuple,
                            names=(self._scanParameter,"modeNr")))

        positionTuple = (self._dh.rho[self._maxRhoInd], self._dh.z[self._yInd])

        return analyticalGRDataFrame, paramDataFrame, positionTuple, self.uc
    #}}}
#}}}
