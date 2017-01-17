#!/usr/bin/env python

from ..collectAndCalcHelpers import (collectSteadyN        ,\
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
from boututils.datafile import DataFile
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
        self._dh = DimensionsHelper(path, uc)
    #}}}

    #{{{_calcPecseliSemiAnalytical
    def _calcPecseliSemiAnalytical(self, path, m, kz, yInd):
        #{{{docstring
        """
        Calculates the analytical growth rates based on the dispersion
        relation assuming Ti = 0 in

        Pecseli, Hans
        Low Frequency Waves and Turbulence in Magnetized Laboratory Plasmas and in the Ionosphere
        IOP Publishing 2016

        Gets the value of the max gradient position and the max gradient
        from the simulations.

        NOTE:
            * Assumes singly ionized particles
            * Assumes relationship between number of nodes and
              wavelength.

        Parameters
        ----------
        path : str
            Path to collect from
        m : int
            Mode number
        kz : float
            The wavenumber of the fluctuations in z.
            Equals 2*pi/lambda_z measured in m^-1.
            Although not a proper wavelength, this can be inferred from
            the 2D poloidal field plots of the fluctuations.
        yInd : int
            The parallel index to slice the data

        Returns
        -------
        om : complex
            The omega in exp(-i[k*x-omega*t])
        """
        #}}}

        # Get the position of the largest gradient
        indMaxGrad = findLargestRadialGradN(path)
        rhoMax = self._dh.rho[indMaxGrad]
        if self._maxRhoInd is None:
            self_maxRhoInd = rhoMax

        ky = m/rhoMax

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
        omStar = calcOmStar(ky, uDE)

        # Calculation of b
        rhoS = self.uc.getNormalizationParameter("rhoS")
        b    = calcPecseliB(ky, rhoS)

        # Calculation of sigmaPar
        omCE = calcOmCE(B)
        with DataFile(os.path.join(path,"BOUT.dmp.0.nc")) as f:
            nuEI = f.read("nuEI")
        sigmaPar = calcSigmaPar(ky, kz, omCE, omCI, nuEI)

        # Test as a function of B
        # For different mode numbers
        om = pecseliAnalytical(omStar, b, sigmaPar)

        return om
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

        multiTuples = []
        fullDict = {"analyticalGR":[], "angularVelocity":[]}

        loopOver = zip(self._steadyStatePaths)

        # Loop over the folders
        for steadyStatePath in self._steadyStatePaths:
            # Set the units converter and the dimension helper
            self._setUcAndDh(steadyStatePath)
            # Obtain the scan value
            scanValue = getScanValue(scanPaths, self._scanParameter)

            # Assume kz is double of the cylinder heigth (observed in
            # the fluctuations)
            kz = 2*(self._dh.z[-1] + 0.5*(self._dh.z[-1] - self._dh.z[-2])

            for modeNr in range(len(1, nModes+1)):
                # Fill the multiIndexTuple and the dict
                multiTuples.append((scanValue, modeNr))

                om = self._calcPecseliSemiAnalytical(steadyStatePath,\
                                                     modeNr         ,\
                                                     kz             ,\
                                                     self._yInd      )

                fullDict["analyticalGR"].append(om.imag)
                fullDict["angularVelocity"].append(om.real)

        # Make the data frame
        analyticalGRDataFrame =\
            pd.DataFrame(fullDict,\
                         index=pd.MultiIndex.from_tuples(\
                            multiTuples,
                            names=(self._scanParameter,"modeNr")))

        positionTuple = (rho[self._maxRhoInd], z[self._yInd])

        return analyticalGRDataFrame, positionTuple, self.uc
    #}}}
#}}}
