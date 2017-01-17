#!/usr/bin/env python

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
    #}}}

    def calcPecseliAnalytical(self, B, kz):
        #{{{docstring
        """
        Calculates the analytical growth rates based on the dispersion
        relation assuming Ti = 0 in

        Pecseli, Hans
        Low Frequency Waves and Turbulence in Magnetized Laboratory Plasmas and in the Ionosphere
        IOP Publishing 2016

        Parameters
        ----------
        kz : float
            T

        Returns
        -------
        """
        #}}}
        # Get the position of the largest gradient
        ind = findLargestRadialGradN(steadyStatePath)
        rhoMax = self._dh.rho[ind]

        # Use the Ellis approach
        ky = m/rhoMax
        kx = ky

        # Calculation of omStar
        # NOTE: For now: Duplication of work as gradient in already
        #       calculated in findLargestRadialGradN
        n    = collectSteadyN(steadyStatePath, yInd)
        dx   = getUniformSpacing(steadyStatePath, "x")
        dndx = DDX(n, dx)
        Te   = self.uc.getParamFromNormDict("Te0")
        uDe  = calcUDE(Te, B, n, dndx)
        omStar = calcOmStar(ky, uDE)

        # Calculation of b
        rhoS = self.uc.getParamFromNormDict("rhoS")
        b = calcPecseliB(ky, rhoS)

        # Calculation of sigmaPar
        sigmaPar = calcSigmaPar(ky, kz, omCE, omCI, nuEI)

        # Test as a function of B
        # For different mode numbers
        pecseliAnalytical(omStar, b, sigmaPar)



#}}}
