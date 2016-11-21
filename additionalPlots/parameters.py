#!/usr/bin/env python

"""
Collection of analytic equations
"""

import numpy as np
from numpy import log, sqrt
import scipy.constants as cst

#{{{
class Parameters(object):
    """Class which contains calculation for plasma parameters."""

    #{{{Constructor
    def __init__(self, mi = cst.m_p):
        # Aliases
        self.mi = mi
        self.me = cst.m_e
        self.q  = cst.e
        self.e0 = cst.epsilon_0
        self.pi = np.pi
    #}}}

    #{{{calcCoulombLog
    def calcCoulombLog(self, Te0, n0):
        coulombLog = log(12.0*self.pi*(self.e0*Te0)**(3/2)/(sqrt(n0)*self.q**3))
        return coulombLog
    #}}}

    #{{{calcVThE
    def calcVThE(self, Te0):
        vThE = sqrt(2.0*Te0/self.me)
        return vThE
    #}}}

    #{{{calcVThI
    def calcVThI(self, Ti0):
        vThI = sqrt(2.0*Ti0/self.mi)
        return vThI
    #}}}

    #{{{calcNuEI
    def calcNuEI(self, ne0, coulombLog, vThE):
        nuEI=1.0/(4.0*self.pi)*(self.q**4.0*ne0/(self.e0**2.0*self.me**2.0))*coulombLog*(1.0/vThE**3.0)
        return nuEI
    #}}}

    #{{{calcNuII
    def calcNuII(self, ne0, coulombLog, vThI):
        nuII=1.0/(2.0*self.pi)*(self.q**4.0*ne0/(self.e0**2.0*self.mi**2.0))*coulombLog*(1.0/vThI**3.0)
        return nuII
    #}}}

    #{{{calcEtaI
    def calcEtaI(self, n, Ti, tauI, omCI):
        etaI = {}

        etaI['\eta_{{i,0}}'] = 0.96*n*Ti*tauI
        etaI['\eta_{{i,1}}'] = 3.0*n*Ti/(10.0*(omCI**2.0)*tauI)
        etaI['\eta_{{i,2}}'] = 4.0*etaI['\eta_{{i,1}}']
        etaI['\eta_{{i,3}}'] = n*Ti/(2.0*omCI)
        etaI['\eta_{{i,4}}'] = 2.0*etaI['\eta_{{i,3}}']

        return etaI
    #}}}

    #{{{calcEtaE
    def calcEtaE(self, n, Te, tauE, omCE):
        etaE = {}

        etaE['\eta_{{e,0}}'] = 0.73*n*Te*tauE
        etaE['\eta_{{e,1}}'] = 0.51*n*Te/(omCE**2.0*tauE)
        etaE['\eta_{{e,2}}'] = 4.0*etaE['\eta_{{e,1}}']
        etaE['\eta_{{e,3}}'] = n*Te/(2.0*omCE)
        etaE['\eta_{{e,4}}'] = 2.0*etaE['\eta_{{e,3}}']

        return etaE
    #}}}

    #{{{calcOmCI
    def calcOmCI(self, B0):
        return self.q*B0/self.mi
    #}}}

    #{{{calcOmCE
    def calcOmCE(self, B0):
        return self.q*B0/self.me
    #}}}

    #{{{CalcCS
    def calcCS(self, Te0, Ti0, N):
        return sqrt((Te0+((N+2.0)/N)*Ti0)/self.mi)
    #}}}

    #{{{
    def calcRhoS(self, cS, omCI):
        return cS/omCI
    #}}}
#}}}
