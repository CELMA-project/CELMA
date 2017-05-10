#!/usr/bin/env python

"""Contains common variables and functions for all coordinate systems"""

from collections import OrderedDict
from IPython.display import display
from sympy import symbols
from sympy import Eq, Derivative
from sympy import sin, cos, atan, sqrt

x, y, z = symbols('x, y, z', real=True)
rho = symbols('rho', positive=True)
theta = symbols('theta', real=True)
cartCoord = [x, y, z]
# Using the Clebsch system B = e^3 x e^1 (the order is important)
cylCoord  = [rho, z, theta]

cartMap = OrderedDict([
            (x, rho*cos(theta)),
            (y, rho*sin(theta)),
            (z, z),
            ])

cylMap = OrderedDict([
            (rho,   sqrt(x**2 + y**2)),
            (z,     z),
            (theta, atan(y/x)),
            ])

#{{{poisson
def poisson(a, b):
    """
    Poisson bracket in cylinder geometry

    NOTE: This is derived for the Clebsch system
    """

    return Derivative(a, theta)*Derivative(b, rho) - Derivative(a, rho)*Derivative(b, theta)
#}}}

#{{{displayVec
def displayVec(v, LHS=None):
    """Displays the vector nicely in a notebook"""

    if LHS == None:
        LHS = 'v'
    if type(v).__name__ == 'ClebschVec':
        if v.covariant:
            binder = '_'
        else:
            binder = '^'
    else:
        binder = '_'

    for coord in cylCoord:
        display(Eq(symbols(LHS + binder + str(coord)),\
                           v.__getattribute__(str(coord))))
#}}}
