#!/usr/bin/env python

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

"""
Makes a parametric plot of the boundary condition of the laplace
invertor for different modes.
"""

def plotLaplaceBCOnCylinder(mode               =  5,\
                            lineThroughPointNr = 45,\
                            extension          ="png"):
    """
    Plots the laplace inversion boundary condition on a cylinder

    Parameters
    ----------
    mode : int
        Mode number
    lineThroughPointNr : int
        What point should the red line go diametrically opposite from
    extension : str
        Extension of the plot to save
    """

    fig = plt.figure(figsize=(18,12))
    ax = fig.gca(projection='3d')
    ax.set_axis_off()

    n = 128


    # NOTE: The parametric curve has the coordinates (x(z), y(z), z)
    #       I. e. for each z, there is a corresponding x,y pair

    zero = np.zeros(n)
    # +1 as we are making this periodically
    theta = np.linspace(0, 2*np.pi, n+1)
    # Until last point (periodic point not counted twice)
    theta = theta[:-1]
    x    = np.sin(theta)
    y    = np.cos(theta)
    # Make the sine
    sine = zero + np.sin(mode*theta)
    # Make x spine (actually not a spine, just a line)
    xSpine = np.linspace(-1.5, 1.5, n)
    # Make y spine (actually not a spine, just a line)
    ySpine = xSpine

    ax.plot(xSpine, zero,   zero, "k--", linewidth = 1.5)
    ax.plot(zero,   ySpine, zero, "k--", linewidth = 1.5)
    ax.plot(x,      y,      zero, "-go", linewidth = 3.0)
    ax.plot(x,      y,      sine, "-bo", linewidth = 3.0)

    # Making vertical bars
    for i in range(len(sine)):
        pointX = (x[i], x[i]  )
        pointY = (y[i], y[i]  )
        pointZ = (0   , sine[i])
        ax.plot(pointX, pointY, pointZ, "k", linewidth=1.0)

    # Making the red line
    diametricallyOpposite = lineThroughPointNr + int(n/2)
    redX = (x[lineThroughPointNr]   , x[diametricallyOpposite])
    redY = (y[lineThroughPointNr]   , y[diametricallyOpposite])
    redZ = (sine[lineThroughPointNr], sine[diametricallyOpposite])
    # NOTE: If we wanted the red line to be behind everything, we would
    #       have set zorder=0
    ax.plot(redX, redY, redZ, "r", linewidth=2.0)

    # Set view
    ax.view_init(elev=55, azim=-52)
    # Zoom in a bit to aviod white space
    ax.dist = 6.5

    fig.tight_layout(pad=0.0)
    fig.savefig("mode_{}.{}".format(mode, extension), transparent=True)
    plt.show()

if __name__ == "__main__":
    # Plot the odd mode
    plotLaplaceBCOnCylinder(mode=5, lineThroughPointNr=57, extension="pdf")
    # Plot the even mode
    plotLaplaceBCOnCylinder(mode=4, lineThroughPointNr=40, extension="pdf")
