#!/usr/bin/env python

"""
Python file which contains:

Analytical expression for growth rates in drift waves given by

   Ellis, R F and Marden-Marshall, E and Majeski, R
   Collisional drift instability of a weakly ionized argon plasma
   Plasma Physics 1980, 22, 2
   doi:10.1088/0032-1028/22/2/002

   Pecseli, Hans
   Low Frequency Waves and Turbulence in Magnetized Laboratory Plasmas and in the Ionosphere
   IOP Publishing 2016
"""

import scipy.constants as cst
import numpy as np

#{{{ellisAnalytical
def ellisAnalytical(omStar, bEllis, nuPar, om1, nuIN):
    #{{{docstring
    """
    Function which calculates the real and imaginary omega from the
    dispersion relation given in the paper by Ellis et al.

    NOTE: Underlying assumption neutral collision dominates Coloumb
          collisions.

    NOTE: All input parameters must be in non-normalized units.

    NOTE: Underlying assumption kx >> (1/n)(dn/dx)

    Parameters
    ----------
    omStar : float
        Parameter telling about electron diamagnetic frequency.
    bEllis : float
        The Ellis definition of b.
        Tells something about the size of the perturbation as compared
        to rhoS.
    nuPar : float
        Parameter which tells something about the conductivity.
    om1 : float
        Parameter which tells something about the streaming of the
        electrons.
    nuIN : float
        The ion-neutral collision frequency

    Returns
    -------
    om : complex
        The omega in exp(-i[k*x-omega*t])
    """
    #}}}
    b = bEllis

    real = omStar/(1+b)
    imag =    (omStar/(nuPar*(1+b)))*\
              ((omStar/((1+b)**2)) + om1)\
            - (b/(b+1))*nuIN

    om = complex(real, imag)
    return om
#}}}

#{{{pecseliAnalytical
def pecseliAnalytical(omStar, bPecseli, sigmaPar):
    #{{{docstring
    """
    Function which calculates the real and imaginary omega from the
    dispersion relation assuming Ti = 0 by Pecseli.

    NOTE: Underlying assumption that sigmaPar/omStar is large and b is
          small.

    NOTE: The rest of the underlying assumption found in section 5.1 in draft.

    NOTE: All input parameters must be in non-normalized units.

    Parameters
    ----------
    omStar : float
        Parameter telling about electron diamagnetic frequency.
    bPecseli : float
        The Pecseli definition of b.
        Tells something about the size of the perturbation as compared
        to rhoS.
    sigmaPar : float
        Parameter which tells something about the conductivity.

    Returns
    -------
    om : complex
        The omega in exp(-i[k*x-omega*t])
    """
    #}}}

    b = bPecseli
    if b >= sigmaPar/omStar:
        message = "Model breaks down as b => sigmaPar/omStar."
        raise RuntimeError(message)

    return complex(omStar*(1-b), omStar**2/sigmaPar)
#}}}

#{{{calcNuPar
def calcNuPar(kz, Te, nuEN):
    #{{{docstring
    """
    Calculates nuPar

    Parameters
    ----------
    kz : float
        The inverse wavelength in z measured in m^-1.
    Te : float
        The electron temperature measured in J.
    nuEN :float
        The electron neutral collision frequency measured in s^-1.

    Returns
    -------
    nuPar : float
        A quantity describing the parallel conductivity
    """
    #}}}
    return kz**2*Te/(cst.m_e*nuEN)
#}}}

#{{{calcOm1
def calcOm1(kz, u0):
    #{{{docstring
    """
    Calculates om1

    Parameters
    ----------
    kz : float
        The inverse wavelength in z measured in m^-1.
        Equals 2*pi/lambda_z.
        NOTE: Typo in the article where it says k2
    u0 : float
        The streaming of electrons in the z direction
    """
    #}}}
    # Calculation of omega1
    return kz*u0
#}}}

#{{{calcSigmaPar
def calcSigmaPar(ky, kz, omCE, omCI, nuEI):
    #{{{docstring
    """
    Calcualtes the parallel sigma

    NOTE: In order to compare with cylindrical geometry, one makes the
          substitution
          x -> rho
          y -> rho*theta (y=rho*sin(theta))

    Parameters
    ----------
    ky : float
        The inverse wavelength in y.
        Equals 2*pi/lambda_y measured in m^-1.
    kz : float
        The inverse wavelength in z.
        Equals 2*pi/lambda_y measured in m^-1.
    omCE : float
        Electron cyclotron frequency measured in s^-1
    omCI : float
        Ion cyclotron frequency measured in s^-1
    nuEI : float
        Electron-ion Coloumb collision frequency measured in s^-1

    Returns
    -------
    sigmaPar : float
        Variable defined as ((kz/ky)**2)*(omCE/nuEI)*omCI.
    """
    #}}}
    return ((kz/ky)**2)*(omCE/nuEI)*omCI
#}}}

#{{{calcEllisB
def calcEllisB(kx, ky, rhoS):
    #{{{docstring
    """
    Calculates the b

    NOTE: In order to compare with cylindrical geometry, one makes the
          substitution
          x -> rho
          y -> rho*theta (y=rho*sin(theta))

    Parameters
    ----------
    kx : float
        The inverse wavelength in x measured in m^-1.
        Equals 2*pi/lambda_x.
    ky : float
        The inverse wavelength in y measured in m^-1.
        Equals 2*pi/lambda_y.
    rhoS : float
        cs/omCI measured in m.

    Returns
    -------
    b : float
       Quantity defined as (kx**2+ky**2)*rhoS**2.
    """
    #}}}

    return (kx**2+ky**2)*rhoS**2
#}}}

#{{{calcPecseliB
def calcPecseliB(ky, rhoS):
    #{{{docstring
    """
    Calculates the b

    NOTE: In order to compare with cylindrical geometry, one makes the
          substitution
          x -> rho
          y -> rho*theta (y=rho*sin(theta))

    Parameters
    ----------
    ky : float
        The inverse wavelength in y measured in m^-1.
        Equals 2*pi/lambda_y.
    rhoS : float
        cs/omCI measured in m.

    Returns
    -------
    b : float
       Quantity defined as (ky**2)*rhoS**2.
    """
    #}}}

    return (ky**2)*rhoS**2
#}}}

#{{{calcOmStar
def calcOmStar(ky, uDE):
    #{{{docstring
    """
    Calculates the omega^*

    NOTE: In order to compare with cylindrical geometry, one makes the
          substitution
          x -> rho
          y -> rho*theta (y=rho*sin(theta))

    Parameters
    ----------
    ky : float
        The inverse wavelength in y measured in m^-1.
        Equals 2*pi/lambda_y.
    uDE : float
        The electron diamagnetic velocity.
        Measured in ms^-1.

    Returns
    -------
    omStar : float
       Quantity defined as ky*uDE.
    """
    #}}}

    return ky*uDE
#}}}

#{{{calcUDE
def calcUDE(Te, B, n, dndx):
    #{{{docstring
    """
    Calculates the electron diamagnetic velocity, assuming constant
    electron temperature.

    Parameters
    ----------
    Te : float
        Constant electron temperature measured in J.
    B : float
        Magnetic field measured in T.
    n : float
        Density measured in m^-3.
    dndx : float
        The density gradient measured in m^-4.

    Returns
    -------
    uDE : float
        Electron diamagnetic velocity (assuming constant Te).
    """
    #}}}

    return -(Te/(cst.e*B))*(1/n)*dndx
#}}}

#{{{calcRhoS
def calcRhoS(cs, omCI):
    #{{{docstring
    """
    Calculates rhoS

    Parameters
    ----------
    cs : float
        The ion sound speed.
    omCI : float
        The ion cyclotron frequency.

    Returns
    -------
    rhoS : float
        rhoS defined as cs/omCI
    """
    #}}}

    return cs/omCI
#}}}

#{{{calcCS
def calcCS(Te, mi):
    #{{{docstring
    """
    Calculates the ion sound speed.

    Parameters
    ----------
    Te : float
        The electron temperature measured in J
    mi : float
        The ion mass in kg,

    Returns
    -------
    cs : float
        The ion sound speed in ms^-1.
    """
    #}}}

    return np.sqrt(Te/mi)
#}}}

#{{{calcOmCI
def calcOmCI(B, mi):
    #{{{docstring
    """
    Calculates the ion cyclotron frequency

    NOTE: Assumes singly ionized particles

    Parameters
    ----------
    B : float
        The magnetic field in T.
    mi : float
        The ion mass in kg.

    Returns
    -------
    omCI : float
        The ion cyclotron frequency in s^-1.
    """
    #}}}

    return cst.e*B/mi
#}}}

#{{{calcOmCE
def calcOmCE(B):
    #{{{docstring
    """
    Calculates the electron cyclotron frequency

    Parameters
    ----------
    B : float
        The magnetic field measured in T.

    Returns
    -------
    omCE : float
        The electron cyclotron frequency in s^-1
    """
    #}}}

    return cst.e*B/cst.m_e
#}}}
