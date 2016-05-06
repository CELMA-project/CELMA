# Generate MMS solutions for the laplace inversion
from boutdata.mms import x, y, z, t
from boutdata.mms import exprToStr
from boutdata.mms import Metric
from boututils.options import BOUTOptions

from sympy import nsimplify
from sympy.utilities.lambdify import lambdify

import numpy as np

import matplotlib.pylab as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import re

#{{{get_metric
def get_metric():
    """Returns the metric"""
    # Create an instance
    metric = Metric()

    # Set the differencing terms (set by a prime by default)
    metric.x = x
    metric.y = y
    metric.z = z

    # Set the metrics
    # Contravariant elements
    metric.g11  = 1
    metric.g22  = 1
    metric.g33  = 1/(x**2)
    metric.g12  = 0
    metric.g23  = 0
    metric.g13  = 0
    # Covariant elements
    metric.g_11 = 1
    metric.g_22 = 1
    metric.g_33 = x**2
    metric.g_12 = 0
    metric.g_23 = 0
    metric.g_13 = 0
    # The Jacobian
    metric.J    = x

    return metric
#}}}

#{{{set_plot_style
def set_plot_style():
    """ Defines the style for the plots"""

    title_size = 30
    plt.rc("font", size = 30)
    plt.rc("axes", labelsize = 25, titlesize = title_size)
    plt.rc("xtick", labelsize = 25)
    plt.rc("ytick", labelsize = 25)
    plt.rc("legend", fontsize = 30)
    plt.rc("lines", linewidth = 3.0)
    plt.rc("lines", markersize = 20.0)
    plt_size = (10, 7)
    n_grid_lines = 400

    return plt_size, n_grid_lines
#}}}

#{{{make_plot
def make_plot(folder, the_vars, plot3d = True, plot2d = False, direction='x',\
              include_aux = False, save = False):
    """3D plot of the MMS variables"""

    plt_size, n_grid_lines = set_plot_style()

    #{{{ Get options from the BOUT.inp file to set x_ax_start and x_ax_end
    myOpts = BOUTOptions(folder)
    # Get variables from the input file
    Lx  = eval(myOpts.geom['Lx'])
    MXG = eval(myOpts.root['MXG'])
    nx  = eval(myOpts.mesh['nx'])

    # Find dx
    dx = Lx/(nx-2*MXG )
    x_ax_start = dx
    x_ax_end   = Lx - dx
    #}}}

    if direction == 'y':
        #{{{ Get options from the BOUT.inp file to set y_ax_start and y_ax_end
        # Get variables from the input file
        Ly  = eval(myOpts.geom['Ly'])
        ny  = eval(myOpts.mesh['ny'])

        # Find dy (remember that ny contains all the gridpoints
        # [including ghosts])
        dy = Ly/ny
        y_ax_start = dy
        y_ax_end   = Ly - dy
        #}}}

    fig_no = 0
    # Levels
    N = 200

    #{{{ Do the plotting
    for cur_var_key in the_vars.keys():
        if not("aux_" in cur_var_key) or\
           (("aux_" in cur_var_key) and include_aux):
            if direction == 'x':
                cur_plt = lambdify((x,z), the_vars[cur_var_key],\
                                   modules = ['numpy'])
                # We would like the z-direction on the y-axis
                y_ax_len = np.linspace(0, 2*np.pi, n_grid_lines)
            elif direction == 'y':
                cur_plt = lambdify((x,y), the_vars[cur_var_key],\
                                   modules = ['numpy'])
                # We would like the y-direction on the y-axis
                y_ax_len = np.linspace(y_ax_start, y_ax_end, n_grid_lines)

            x_ax_len = np.linspace(x_ax_start, x_ax_end, n_grid_lines)
            X_ax_len, Y_ax_len = np.meshgrid(x_ax_len, y_ax_len)

            if plot3d:
                # Plot the variables in a 3D plot
                fig_no += 1
                fig = plt.figure(fig_no, figsize = plt_size)
                ax1 = fig.add_subplot(111)
                # Plot the plot
                ax1 = fig.gca(projection='3d')
                ax1.plot_surface(X_ax_len, Y_ax_len, cur_plt(X_ax_len, Y_ax_len),\
                                 cmap = cm.RdYlBu_r,\
                                 linewidth = 0)
                # Set the labels
                ax1.set_xlabel('x')
                if direction == 'x':
                    ax1.set_ylabel('z')
                elif direction == 'y':
                    ax1.set_ylabel('y')
                ax1.set_zlabel(cur_var_key)
                # Set the grid
                ax1.grid()

            if plot2d:
                if direction == 'x':
                    # Transform to cylindrical
                    X_RT = X_ax_len*np.cos(Y_ax_len)
                    Y_RT = X_ax_len*np.sin(Y_ax_len)
                elif direction == 'y':
                    X_RT = X_ax_len
                    Y_RT = Y_ax_len
                # Plot the variables in a 2D plot
                fig_no += 1
                fig = plt.figure(fig_no, figsize = plt_size)
                ax2 = fig.add_subplot(111)
                # Plot the plot
                # zorder decides what should be drawn first
                cont = ax2.contourf(X_RT, Y_RT, cur_plt(X_ax_len, Y_ax_len),\
                                    N, cmap = cm.RdYlBu_r, zorder=-20)
                cbar = plt.colorbar(cont)
                cbar.ax.set_ylabel(cur_var_key)
                # Set the labels
                ax2.set_xlabel(r'$\rho$')
                if direction == 'x':
                    ax2.set_ylabel(r'$\rho$')
                elif direction == 'y':
                    ax2.set_ylabel(r'$z$')
                # Set the grid
                ax2.grid()

                if save:
                    # Set zorder value below which artists will be rasterized
                    # If this is not set, everything will be vecotrized, giving
                    # plots with a huge size
                    ax2.set_rasterization_zorder(-10)
                    plt.tight_layout()
                    plt.savefig('{}.pdf'.format(cur_var_key))
    #}}}
    # Show the plots
    plt.show()
#}}}

#{{{Cast to string, replace and print
def BOUT_print(the_vars, rational=False):
    """Print the variables in BOUT++ style"""

    for cur_var_key in the_vars.keys():
        if rational:
            the_vars[cur_var_key] = nsimplify(the_vars[cur_var_key])
        the_vars[cur_var_key] = exprToStr(the_vars[cur_var_key])
        the_vars[cur_var_key] = re.sub(r'\b' + "x" + r'\b', "geom:xl", the_vars[cur_var_key])
        print("\n[{}]".format(cur_var_key))
        print(the_vars[cur_var_key])
#}}}
