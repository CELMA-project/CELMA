#!/usr/bin/env python

from ..collectAndCalcHelpers import (collectSteadyN        ,\
                                     DDX                   ,\
                                     getUniformSpacing     ,\
                                     findLargestRadialGradN,\
                                    )
from ..superClasses import CollectAndCalcSuperClass
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
class CollectAndCalcAnalyticGrowthRates(CollectAndCalcSuperClass):
    """
    Class for collecting and calculating variables needed for
    calculations of the analytical growth rates.
    """

    #{{{constructor
    def __init__(self, *args, **kwargs):
        #{{{docstring
        """
        This constructor will:
            * Call the parent constructor

        NOTE: Mode is set to normal as the poloidal average equals the
              zeroth mode.

        Parameters
        ----------
        *args : positional arguments
            See parent constructor for details.
        *kwargs : keyword arguments
            See parent constructor for details.
        """
        #}}}

        # Call the constructor of the parent class
        super().__init__(*args, **kwargs)

        self._path = self._collectPaths[0]
    #}}}

    #{{{calcPecseliSemiAnalytical
    def calcPecseliSemiAnalytical(self, m, kz, yInd):
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
        indMaxGrad = findLargestRadialGradN(self._path)
        rhoMax = self._dh.rho[indMaxGrad]

        ky = m/rhoMax

        # Calculate B
        mi   = self.uc.getNormalizationParameter("mi")
        omCI = self.uc.getNormalizationParameter("omCI")*mi/cst.e
        B    = omCI*mi/cst.e

        # Calculation of omStar
        # NOTE: For now: Duplication of work as gradient in already
        #       calculated in findLargestRadialGradN
        n    = collectSteadyN(self._path, yInd)
        dx   = getUniformSpacing(self._path, "x")
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
        with DataFile(os.path.join(self._path,"BOUT.dmp.0.nc")) as f:
            nuEI = f.read("nuEI")
        sigmaPar = calcSigmaPar(ky, kz, omCE, omCI, nuEI)

        # Test as a function of B
        # For different mode numbers
        om = pecseliAnalytical(omStar, b, sigmaPar)

        return om
    #}}}
#}}}
