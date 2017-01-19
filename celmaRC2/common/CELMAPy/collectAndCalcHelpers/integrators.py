#!/usr/bin/env python

"""
Contains functions for taking poloidal averages
"""

import numpy as np

#{{{parallelIntegration
def parallelIntegration(f, dy, startInd=None, endInd=None):
    #{{{docstring
    """
    Returns the parallel integration of a field.

    The integration is approximated by the

    Parameters
    ----------
    f : array-4d
        The field to find the poloidal average of.
        The field must be a 4D field.
    dy : float
        The gridspacing in y.
    startInd : [None|int]
        If not None, the parallel index to start the integration from.
    endInd : [None|int]
        If not None, the parallel index to end the integration at.

    Returns
    -------
    out : array-4d
        The parallel integrated of the field.
    """
    #}}}

    out = f[:,:,startInd:endInd,:].sum(axis=2)
    # Expand to a 4d numpy array
    out = np.expand_dims(out*dy, axis=2)

    return out
#}}}

#{{{poloidalIntegration
def poloidalIntegration(f, rho, startInd=None, endInd=None):
    #{{{docstring
    """
    Returns the parallel integration of a field.

    Parameters
    ----------
    f : array-4d
        The field to find the poloidal average of.
        The field must be a 4D field, and should not include the last
        poloidal slice in order not to count the same point twice (i.e.
        the domain should go from [0,2pi[)
    rho : array-1d
        In a cylinder geometry, the poloidal line element is spanned by
        rho*dtheta.
        Must have the same dimension at axis 1 of f.
    startInd : [None|int]
        If not None, the poloidal index to start the integration from.
    endInd : [None|int]
        If not None, the poloidal index to end the integration at.

    Returns
    -------
    out : array
        The poloidal integrated of the field.
    """
    #}}}

    # Obtain dTheta
    # NOTE: Theta in [0, 2*pi] so we must include the last point
    dTheta = 2*np.pi/f.shape[-1]
    out = f[:,:,:,startInd:endInd].sum(axis=-1)
    # Expand to a 4d numpy array
    out = np.expand_dims(out*rho*dTheta, axis=-1)

    return out
#}}}

#{{{radialIntegration
def radialIntegration(f, dx, startInd=None, endInd=None):
    #{{{docstring
    """
    Returns the parallel integration of a field.

    The integration is approximated by the

    Parameters
    ----------
    f : array-4d
        The field to find the poloidal average of.
        The field must be a 4D field.
    dy : float
        The gridspacing in y.
    startInd : [None|int]
        If not None, the parallel index to start the integration from.
    endInd : [None|int]
        If not None, the parallel index to end the integration at.

    Returns
    -------
    out : array-4d
        The parallel integrated of the field.
    """
    #}}}

    out = f[:,startInd:endInd,:,:].sum(axis=1)
    # Expand to a 4d numpy array
    out = np.expand_dims(out*dx, axis=1)

    return out
#}}}
