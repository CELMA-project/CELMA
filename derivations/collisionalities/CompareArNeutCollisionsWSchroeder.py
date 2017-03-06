#!/usr/bin/env python

"""
Contains functions for calculation of the ion-neutral collision
frequency for Argon.

NOTE: Our calculation does not equals the ones found in table 5.1 in
Shroeder, C. - Ph.D. Thesis, 2003, but are of the same order of
magnitudes.
"""

import numpy as np
import scipy.constants as cst

#{{{crossSecArCELowEn
def crossSecArCELowEn(TiEV):
    #{{{docstring
    """
    Calculates the cross section for charge exchange in Argon for low
    energies.

    Source:
    A. Anders. A Formulary for Plasma Physics. Akademie-Verlag, Berlin, 1990.
    Equation (B.7) from Shroeder, C. - Ph.D. Thesis, 2003

    Parameters
    ----------
    TiEV : float
        The ion temperature measured in electron volts.

    Returns
    -------
    crossSecArCE : float
        Cross section for charge exchange in Ar measured in m^2
    """
    #}}}
    return 4.8e-19*(1+0.14*np.log(1/TiEV))**2.0
#}}}

#{{{nuINCELowEn
def nuINCELowEn(nn, Ti):
    #{{{docstring
    """
    Calculated the charge exchange collision frequency for Ar at low
    energies.

    Source:
    Loeiten, M - Ph.D. Thesis, 2017

    Parameters
    ----------
    nn : float
        The neutral density.
    Ti : float
        The Ar temperature measured in Joules.

    Returns
    -------
    nuINCE : float
        The charge exchange collision frequency for Ar measured in s^-1
    """
    #}}}
    mi = 39.948*cst.m_u
    nuINCE = (8*np.sqrt(2)/3)*(nn/np.sqrt(np.pi))*(np.sqrt(Ti/mi))*crossSecArCELowEn(Ti/cst.eV)
    return nuINCE
#}}}

#{{{nuINTotLowEn
def nuINTotLowEn(nn, Ti):
    #{{{docstring
    """
    Calculates the total ion-neutral collision frequencies for Ar at low
    energies.

    We have that the charge exchange and elastic collisions are almost
    equal and dominating over other processes at low energies.

    Sources:
    Lieberman, M.A. and Lichtenberg, A.J. - Principles of Plasma Discharges and Materials Processing, Wiley, 2005
    Figure B.4 and B.5 in Shroeder, C. - Ph.D. Thesis, 2003

    Parameters
    ----------
    nn : float
        The neutral density.
    Ti : float
        The Ar temperature measured in Joules.

    Returns
    -------
    nuINTot : float
        The total ion-neutral collision frequencies for Ar measured in
        s^-1
    """
    #}}}

    return 2.0*nuINCELowEn(nn, Ti)
#}}}

if __name__ == "__main__":
    # Reference density used in page 65 in Shroeder, C. - Ph.D. Thesis, 2003
    nn = 2.5e19

    # Reference temperature used in page 65 in Shroeder, C. - Ph.D. Thesis, 2003
    TiRoom = 273*cst.Boltzmann

    # Collision frequency used in page 65 in Shroeder, C. - Ph.D. Thesis, 2003
    nuINS = 6e4

    nuINTot = nuINTotLowEn(nn, TiRoom)
    message = ("Our frequency: {}, Shroeder, C. frequency {}\nSchroeder is a "
               "factor {} larger").\
                format(nuINTot, nuINS, nuINS/nuINTot)
    print(message)
