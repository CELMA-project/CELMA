#!/usr/bin/env python

"""
Contains functions for estimation of the electron-neutral collision
frequency for Argon.

NOTE:
The collision frequency is found from a polynomial fit of figure B.3 in [1]
Figure B.3 is obtained from numerical integration of the cross section
found in [2], which we have not succeeded to obtain.
[1] Also gives a plot of the cross section in figure B.2, but numerical
integration of this gives poor agreement with what is found in B.3.
The data from figure B.3 has been extracted using [3].
The data should only be used in the interpolation range.
Extrapolation would be more appropriate if a "natural spline" was used,
but it is trickier to implement as it yields one equation per segment.

[1] Shroeder, C. - Ph.D. Thesis, 2003
[2] M. Hayashi. - Recommended values of transport cross sections for elastic collision and total collision cross section for electrons in atomic and molecular gases.  Institute of Plasma Physics, Nagoya University, Report IPPJ-AM-19, 1981
[3] Rohatgi, Ankit - WebPlotDigitalizer available at http://arohatgi.info/WebPlotDigitizer/
"""

import numpy as np
import matplotlib.pylab as plt
import scipy.constants as cst
from scipy.integrate import simps

# The data collected using [3]
#{{{nuENTe
nuENTe =\
(
(0.15318610095220686, 40010.034955955125),\
(0.12169581045844202, 31979.080924571746),\
(0.25395503053225266, 58684.79503080276),\
(0.3200846405691582, 81829.63590014179),\
(0.405108424902322, 117449.41697193087),\
(0.5090263835317446, 166155.03331402616),\
(0.6318385164574258, 234573.7101938146),\
(0.7735448236793658, 321894.942201736),\
(0.9435923923456939, 440810.7400617494),\
(1.14198122245641, 594407.872866794),\
(1.3498171397152556, 769024.448291626),\
(1.5576530569741007, 960256.0717290965),\
(1.7654889742329458, 1155725.2625334312),\
(1.973324891491791, 1354903.344757911),\
(2.181160808750636, 1555360.1026345377),\
(2.388996726009481, 1759851.3438690784),\
(2.596832643268326, 1957501.1252043194),\
(2.804668560527172, 2143284.163094236),\
(3.0125044777860173, 2331331.493256673),\
(3.2203403950448624, 2506063.54949089),\
(3.4281763123037075, 2672735.476074832),\
(3.6360122295625525, 2835548.692500537),\
(3.8438481468213985, 2988579.2836881336),\
(4.051684064080243, 3129240.877690575),\
(4.259519981339089, 3272220.1619228777),\
(4.467355898597933, 3394860.13749368),\
(4.675191815856779, 3512852.1683284785),\
(4.883027733115625, 3630171.717188889),\
(5.0908636503744695, 3731742.741574129),\
(5.2986995676333155, 3826086.9999922463),\
(5.5065354848921615, 3912520.28082794),\
(5.714371402151006, 3995652.1304217908),\
(5.922207319409852, 4096668.4285019743),\
(6.130043236668696, 4145405.709159108),\
(6.337879153927542, 4233485.826112874),\
(6.545715071186386, 4283850.797178092),\
(6.753550988445232, 4357659.877037704),\
(6.961386905704078, 4397928.580554107),\
(7.169222822962922, 4485476.126272952),\
(7.377058740221768, 4485476.126272952),\
(7.584894657480612, 4532878.558794497),\
(7.792730574739458, 4617041.404581049),\
(8.000566491998304, 4617041.404581049),\
(8.208402409257149, 4629191.5620582625),\
(8.416238326515995, 4727550.990011621),\
(8.624074243774839, 4752465.676219845),\
(8.831910161033685, 4752465.676219845),\
(9.03974607829253, 4752465.676219845),\
(9.247581995551375, 4752465.676219845),\
(9.455417912810221, 4802689.649592133),\
(9.663253830069065, 4891862.130852429),\
(9.871089747327911, 4927345.477252272),\
)
#}}}

TeEV = np.array([val[0] for val in nuENTe])
nuEN = np.array([val[1] for val in nuENTe])

#{{{getPolyFit
def getPolyFit(deg = 6):
    #{{{docstring
    """
    Get the polynomial fit

    Parameters
    ----------
    deg : int
        Degree of polynomial to use

    Returns
    -------
    p : array
        The polynomial coefficients in ascending order.
        NOTE: This is revesered as compared to np.polyfit
    polyStr : str
        String of the polynomial function.
        Meant to use in a C++ program.
    """
    #}}}
    cmap = plt.get_cmap("viridis")
    p = np.polyfit(TeEV, nuEN, deg)
    # polyfit gives coefficients of highest polynomial first
    p = p[::-1]
    polyStr = ""
    for i in range(deg+1):
        polyStr += "{:+f}*pow(Te0_, {})\n".format(p[i], i)

    if polyStr[0] == "+":
        polyStr = polyStr[1:]
    polyStr = "nuEN = (nn/2.5e-19)*(\n" + polyStr
    polyStr += ");"

    return p, polyStr
#}}}

#{{{evalPolynomial
def evalPolynomial(p, x):
    #{{{docstring
    """
    Evaluate the polynomial in x

    Parameters
    ----------
    p : array
        The polynomial coefficients in ascending order.
        NOTE: This is revesered as compared to np.polyfit.
    x : array-like
        Value to evaluate to polynomial in.

    Returns
    -------
    y : array-like
        The polynomial evaluated in x
    """
    #}}}
    deg = len(p)
    y = np.zeros(len(x))
    for i in range(deg):
        y += p[i]*x**i

    return y
#}}}

if __name__ == "__main__":
    # Get the polynomial and y
    p, polyStr = getPolyFit()
    y = evalPolynomial(p, TeEV)

    # Plot
    fig, ax = plt.subplots()
    ax.plot(TeEV, nuEN, lw=7, c="gray", label="Schroeder")
    ax.plot(TeEV,    y, lw=4, ls="--", c="k", label="Interpolation")
    ax.set_yscale("log")
    ax.set_xlabel(r"$T_e$ $[\mathrm{eV}]$")
    ax.set_ylabel(r"$\nu_{en}$ $[\mathrm{s}^{-1}]$")
    ax.grid(b = True, which="major")
    ax.grid(b = True, which="minor", axis="y")
    ax.legend(fancybox = True, framealpha = 0.5, loc = "best")
    plt.show()

    print("Integrating polynomial:")
    print(polyStr)
