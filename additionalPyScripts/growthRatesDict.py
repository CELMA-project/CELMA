#!/usr/bin/env python

"""
Python file which contains:

1. Analytical expression for growth rates in drift waves given by

   Ellis, R F and Marden-Marshall, E and Majeski, R
   Collisional drift instability of a weakly ionized argon plasma
   Plasma Physics 1980, 22, 2
   doi:10.1088/0032-1028/22/2/002

2. Data of growth rates obtained from
    Figure 7.2 at page 93 of
    Ph.D. thesis
    Schr{\"{o}}der, C.
    Experimental investigations on drift waves in linear magnetized plasmas
    Greifswald 2002

    These growth rates are obtained by pixel counting using

    WebPlotDigitizer - http://arohatgi.info/WebPlotDigitizer/app/
"""

import scipy.constants as cst
import numpy as np

#{{{ellisAnalytical
def ellisAnalytical(m_i, Te, Ti, B, n, dndx,\
                    k_x, k_y, k_z, nu_en, nu_in, k_2, u_0):
    #{{{docstring
    """
    Function which calculates the real and imaginary omega from the
    dispersion relation given in the paper by Ellis et al.

    NOTE: All input parameters must be in non-normalized units.

    NOTE: Underlying assumption k_x >> (1/n)(dn/dx)

    NOTE: In order to compare with cylindrical geometry, one makes the
          substitution
          x -> rho
          y -> rho*theta (y=rho*sin(theta)

    Parameters
    ----------
    m_i : float
        Ion mass
    Te : float
        Electron temperature in eV
    Ti : float
        Ion temperature in eV
    B : float
        Magnetic field in T
    n : float
        Plasma density
    dndx : float
        The x derivative of n at the given point
    k_x : float
        The inverse wavelength in x (2*pi/lambda)
    k_y : float
        The inverse wavelength in y (2*pi/lambda)
    k_z : float
        The inverse wavelength in z (2*pi/lambda)
    nu_en : float
        Electron neutral collision rates
    nu_in : float
        Ion neutral collision rates
    k_2 : float
        Not specified in the article, but suspect it is a typo and
        should have been k_z
    u_0 : float
        Electron streaming in the z direction


    Returns
    -------
    imag : float
        Im(omega) in exp(-i[k*x-omega*t])
    real : float
        Re(omega) in exp(-i[k*x-omega*t])
    """
    #}}}

    # Conversion from eV to J
    TeJ = Te*cst.e
    TiJ = Ti*cst.e
    # Calculation of the electron diamagnetic drift
    v_d = -(Te/(cst.e*B))*(1/n)*(dndx)
    # Calculation of om_star (drift wave frequency)
    om_star = k_y*v_d

    # Calculation of c_s
    N = 3 # Degrees of freedom
    c_s = np.sqrt((Te0J+((N+2.0)/N)*Ti0J)/mi)
    # Calculation of om_ci
    om_ci = cst.e*B/ m_i
    # Calculation of rho_s
    rho_s = c_s/om_ci
    # Calculation of k_perp
    k_perp = np.sqrt(k_x**2+k_y**2)
    # Calculation of b
    b = k_perp**2/(rho_s**2)

    # Calculation of nu_par
    nu_par = k_z**2*TeJ/(cst.m_e*nu_en)

    # Calculation of omega_1
    omega_1 = k_2*u_0

    real = om_star/(1+b)
    imag =    (om_star/(nu_par*(1+b)))*((nu_star/((1+b)**2)) + om_1)\
            - (b/(b+1))*nu_in

    return imag, real
#}}}

#{{{pixelGrowthRates
pixelGrowthRates =\
{\
"Ellis" :\
    {\
    "gaussianRealX" : np.array(\
        [0.9688657696057688,\
         1.9053919021362804,\
         2.937950067667301 ,\
         3.9024784265732135,\
         4.954128009743912 ,\
         5.866409638739788 ,\
         6.92366707356243  ,\
         7.811137332462596 ,\
         8.839309001861743 ,\
         9.803173339404868 ,\
         10.85052351287121 ,\
         11.843356724364352,\
         12.944412805989568,\
         13.839437504389698,\
         14.846194972299564,\
        ]),\
    "gaussianRealY" : 1.0e4*np.array(\
        [7.191581467807193 ,\
         7.323717956958561 ,\
         7.392962464937021 ,\
         7.408592004483947 ,\
         7.333522307038992 ,\
         7.109093052623821 ,\
         6.770324100617999 ,\
         6.384129872500598 ,\
         5.904585133440005 ,\
         5.4977210313879175,\
         5.132227087259571 ,\
         4.81086474539114  ,\
         4.420342226667317 ,\
         4.239098200580987 ,\
         4.036043841899711 ,\
        ]),\
    "gaussianImagX" : np.array(\
        [-0.37343680839158866,\
          0.847737231255755  ,\
          1.9187074479585657 ,\
          2.9692037067195196 ,\
          3.991809500478354  ,\
          5.040281262539673  ,\
          6.036961641570415  ,\
          7.039766091944421  ,\
          8.10223474482926   ,\
          8.997903243788041  ,\
          10.036772630935246 ,\
          11.11776194820958  ,\
          12.08295877528333  ,\
          13.09174291468203  ,\
          14.129849987379963 ,\
          15.12663492169543  ,\
        ]),\
    "gaussianImagY" : 1.0e4*np.array(\
        [-0.3325208996233364 ,\
          0.14827395418368994,\
          0.7529820414958159 ,\
          1.4192093323812642 ,\
          2.1098668771244142 ,\
          2.9178899224861645 ,\
          3.6280134292085204 ,\
          4.2166051265669395 ,\
          4.7048770659517665 ,\
          5.064705436705774  ,\
          5.2663325023592815 ,\
          5.375267613308932  ,\
          5.555290494921675  ,\
          5.574215637762114  ,\
          5.599073662298814  ,\
          5.593719673623449  ,\
        ]),\
    "fullRealX"     : np.array(\
        [0.912256503523649 ,\
         1.9185349692057099,\
         3.890453805070549 ,\
         5.881805058226794 ,\
         7.897726996430939 ,\
         9.852985146542943 ,\
         11.842958488959638,\
         13.950125606332303,\
        ]),\
    "fullRealY"     : 1.0e4*np.array(\
        [6.248043040154126 ,\
         6.321782958342652 ,\
         6.51923991218936  ,\
         6.326086708996286 ,\
         5.761390386537644 ,\
         5.099242391130147 ,\
         4.49271514628842  ,\
         3.8802172383442635,\
        ]),\
    "fullImagX"     : np.array(\
        [0.8731070062504962,\
         1.9167470648442477,\
         3.9837821644793454,\
         6.001239240486946 ,\
         8.084368658081328 ,\
         10.089703051668241,\
         12.062171072777737,\
         14.1810025346796  ,\
        ]),\
    "fullImagY"     : 1.0e4*np.array(\
        [-0.1556811287765072 ,\
          0.38429346573895984,\
          1.5800177886904052 ,\
          2.8919047365519983 ,\
          3.836585250002079  ,\
          4.419589273336161  ,\
          4.6847287120319905 ,\
          4.805646343296656  ,\
        ]),\
    }
"Naulin":\
    {\
    "realX" : np.array(\
        [1 ,\
         2 ,\
         3 ,\
         4 ,\
         5 ,\
         6 ,\
         7 ,\
         8 ,\
         9 ,\
         10,\
         11,\
         12,\
         13,\
         14,\
         15,\
        ]),\
    "realY" : 1.0e4*np.array(\
        [9.001530368698035 ,\
         8.841546067703412 ,\
         8.619276367381557 ,\
         8.307462874371375 ,\
         7.848604838326207 ,\
         7.301575534905215 ,\
         6.683401044951192 ,\
         6.072137554513766 ,\
         5.4823497331513   ,\
         4.945508625708371 ,\
         4.42645647466204  ,\
         3.9830734165289154,\
         3.6057231458194314,\
         3.2370978381929074,\
         2.917762249568579 ,\
        ]),\
    "imagX" : np.array(\
        [1 ,\
         2 ,\
         3 ,\
         4 ,\
         5 ,\
         6 ,\
         7 ,\
         8 ,\
         9 ,\
         10,\
         11,\
         12,\
         13,\
         14,\
         15,\
        ]),\
    "imagY" : 1.0e4*np.array(\
        [0.6324091614809895,\
         1.5259786830596767,\
         2.4407017973107195,\
         3.3999318999704755,\
         4.3066344007095765,\
         5.111342732768572 ,\
         5.7888304197722835,\
         6.263442954574254 ,\
         6.667167426165495 ,\
         6.942132305716113 ,\
         7.126730460970503 ,\
         7.237862036272428 ,\
         7.317362041369735 ,\
         7.42085478450281  ,\
         7.376970364253589 ,\
        ]),\
    }\
}
#}}}
