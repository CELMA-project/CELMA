#!/usr/bin/env python

def get_coeff(f, x):
    """
    Solves for the coefficients in a Netwon polynomial using divided
    differences

    Input:
    f - the value of the function evaluated in x
    x - the corresponding x values

    Output:
    a - a list of the coefficients ordered like this [a0,a1...aN]

    Sources:
    Laney, C.B. - Computational Gasdynamics page 135
    www.math-cs.gordon.edu/courses/ma342/python/divided_difference.py
    """

    # Find the looping index
    n = len(x) - 1
    # Initially copy all the f values to a
    a = f.copy()

    for i in range( 1, n + 1 ):
        for j in range( n, i - 1 , -1 ):
            # Notice that a[j] uses a[j]. This is no problem as
            # initially a[j] is equal to f
            a[j] = (a[j] - a[j-1]) / (x[j] - x[j-i])

    return a

def get_polynomial(a, x, xeval):
    """
    Returns a Netwon polynomial

    Input:
    a     - The coefficients of the polynomial
    x     - The postition where the points are evaluated
    xeval - The position to evaluate the polynomial

    Output:
    p - The interpolating polynomial

    Source:
    Laney, C.B. - Computational Gasdynamics page 135
    """

    # Find the order
    order = len(a) - 1
    # Find the integration polynomial
    p = a[order]
    # Looping backwards
    for i in range( order, 0, -1 ):
        # Notice that we will use the previous value of p to compute p
        p = p * (xeval - x[i-1] ) + a[i-1]

    return p
