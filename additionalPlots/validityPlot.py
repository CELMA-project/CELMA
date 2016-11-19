#!/usr/bin/env python

"""
Scripts which plots the valid region of drift approximation.
"""


# FIXME: y-axis will scale with m_H/m_alpha.  Have this on the y-axis



from scipy.optimize import fsolve
import scipy.constants as cst
import numpy as np
import matplotlib.pylab as plt

me = cst.m_e
q  = cst.e
e0 = cst.epsilon_0
pi = np.pi

# vThE = (2.0*Te0/me)**0.5
# nu_ei/om_ci = (mi*ne*q**3*coulombLog)/(4*B*e0**2*me**2*pi*vThE**3.0)
#
# nu_ei/om_ci = (mi*ne*q**3*log(12.0*pi*(e0*Te0)**(3/2)/(n0*self.q**3)))/(4*B*e0**2*me**2*pi*(2.0*Te0/me)**0.5)

mis = [cst.m_p, 39.948*cst.u]
Bs = [1e-2, 5e-2, 1e-1]
Bs = Bs[::-1]
TeEV = np.linspace(1, 15, 20)
TeJ  = TeEV*cst.eV
# TeJ  = TeEV*cst.eV

# n = np.logspace(1e16, 1e20, 200)

def nuEIOverOmCIMinusOne(n, Te, B, mi):
    shouldBeZero =\
            (mi*n*(q**3.0)*np.log(12.0*pi*((e0*Te)**(3.0/2.0))/(np.sqrt(n)*q**3.0)))/\
            (4.0*B*(e0**2.0)*(me**2.0)*pi*((2.0*Te/me)**(3.0/2.0)))\
            - 1
    return shouldBeZero

seqCMap = plt.get_cmap("viridis")
colors = seqCMap(np.linspace(0, 1, len(Bs)))[::-1]

for mi in mis:
    fig, ax = plt.subplots()

    for BNr in range(len(Bs)):
        n = np.zeros(TeJ.shape)
        for i in range(len(TeJ)):
            n[i] = fsolve(nuEIOverOmCIMinusOne, x0=1e16, args=(TeJ[i], Bs[BNr], mi))

        n = tuple(n)

        ax.plot(TeEV, n, color=colors[BNr], label="B = {}".format(Bs[BNr]))
        ax.fill_between(TeEV, n, facecolor=colors[BNr], edgecolor=colors[BNr])
    ax.grid()
plt.show()
