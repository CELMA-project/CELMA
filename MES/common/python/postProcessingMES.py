#!/usr/bin/env python

"""Post processing which performs MES"""

from boutdata import collect
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
from pylab import plot
import numpy as np
import os

#{{{perform_MES_test
def perform_MES_test(\
                     paths,\
                     extension     = 'png' ,\
                     show_plot     = False ,\
                     xz_error_plot = False ,\
                     xy_error_plot = False ,\
                     use_dx        = True  ,\
                     use_dy        = True  ,\
                     use_dz        = True  ,\
                     y_plane       = None  ,\
                     z_plane       = None  ,\
                     yguards       = False ,\
                     ):
    """Collects the data members belonging to a convergence plot"""

    # Figure out the directions
    directions = {'dx':use_dx, 'dy':use_dy, 'dz':use_dz}

    # Make a variable to store the errors and the spacing
    data = {'error_2':[], 'error_inf':[], 'spacing':[],\
            'error_field':[]}

    # Loop over the runs in order to collect
    for path in paths:
        # Collect n_solution - n_numerical
        error_field = collect('e', path=path, info=False,\
                              xguards = False, yguards = yguards)

        if len(error_field.shape) == 4:
            # Pick the last time point
            error_field = error_field[-1]

        # Add the error field
        Lx = collect('Lx', path=path, info=False)
        dx = collect('dx', path=path, info=False)
        dy = collect('dy', path=path, info=False)
        dz = collect('dz', path=path, info=False)
        data['error_field'].append({'field':abs(error_field),\
                                    'Lx':Lx                 ,\
                                    'dx':dx[0,0]            ,\
                                    'dy':dy[0,0]            ,\
                                    'dz':dz})

        # The error in the 2-norm and infintiy-norm
        data['error_2']  .append( np.sqrt(np.mean( error_field**2.0 )) )
        data['error_inf'].append( np.max(np.abs( error_field )) )

        # We are interested in the max of the spacings, so we make
        # an appendable list
        max_spacings = []
        # Collect the spacings
        if use_dx:
            dx_spacing = collect("dx", path=path, info=False,\
                                 xguards = False, yguards = False)
            max_spacings.append(np.max(dx_spacing))
        if use_dy:
            dy_spacing = collect("dy", path=path, info=False,\
                                 xguards = False, yguards = False)
            max_spacings.append(np.max(dy_spacing))
        if use_dz:
            dz_spacing = collect("dz", path=path, info=False,\
                                 xguards = False, yguards = False)
            max_spacings.append(np.max(dz_spacing))

        # Store the spacing in the data
        data['spacing'].append(np.max(max_spacings))

    # Sort the data
    data = sort_data(data)

    # Find the order of convergence in the 2 norm and infinity norm
    order_2, order_inf = get_order(data)

    # Get the root name of the path (used for saving files)
    root_folder = os.path.join(paths[0].split('/')[0], "MES_results")
    # Create folder if it doesn't exist
    if not os.path.exists(root_folder):
        os.makedirs(root_folder)
        print(root_folder + " created\n")


    # Get the name of the plot based on the first folder name
    name = paths[0].split('/')
    # Remove the root folder and put 'MES' in front
    name = 'MES'

    # Print the convergence rate
    print_convergence_rate(data, order_2, order_inf, root_folder, name)

    # Plot
    # We want to show the lines of the last orders, so we send in
    # order_2[-1] and order_inf[-1]
    do_plot(data, order_2[-1], order_inf[-1],\
            root_folder, name, extension,\
            xz_error_plot, xy_error_plot, show_plot,\
            y_plane, z_plane, directions)
#}}}

# Help functions
#{{{sort_data
def sort_data(data):
    """Sorts the data after highest grid spacing"""

    # Sort the data in case it is unsorted
    list_of_tuples_to_be_sorted =\
        list(zip(data['spacing'], data['error_inf'], data['error_2']))

    # Sort the list
    # Note that we are sorting in reverse order, as we want the
    # highest grid spacing first
    sorted_list = sorted(list_of_tuples_to_be_sorted, reverse = True)
    # Unzip the sorted list
    data['spacing'], data['error_inf'], data['error_2'] =\
            list(zip(*sorted_list))

    return data
#}}}

#{{{get_order
def get_order(data):
    # Initialize the orders
    order_2 = [np.nan]
    order_inf = [np.nan]

    # The order will be found by finding a linear fit between two
    # nearby points in the error-spacing plot. Hence, we must let
    # the index in the for loop run to the length minus one
    for index in range(len(data['spacing']) - 1):
        # p = polyfit(x,y,n) finds the coefficients of a polynomial p(x)
        # of degree that fits the data, p(x(i)) to y(i), in a least squares
        # sense.
        # The result p is a row vector of length n+1 containing the
        # polynomial coefficients in descending powers
        spacing_start   = np.log(data['spacing'][index])
        spacing_end     = np.log(data['spacing'][index + 1])
        error_start_2   = np.log(data['error_2'][index])
        error_end_2     = np.log(data['error_2'][index + 1])
        error_start_inf = np.log(data['error_inf'][index])
        error_end_inf   = np.log(data['error_inf'][index + 1])
        # Finding the order in the two norm
        order = np.polyfit([spacing_start, spacing_end],\
                           [error_start_2, error_end_2], 1)
        # Append it to the order_2
        order_2.append(order[0])

        # Finding the infinity order
        order = np.polyfit([spacing_start, spacing_end],\
                           [error_start_inf, error_end_inf], 1)
        # Append it to the order_inf
        order_inf.append(order[0])

    return order_2, order_inf
#}}}

#{{{print_convergence_rate
def print_convergence_rate(data, order_2, order_inf, root_folder, name):
    "Prints the convergence rates to the screen and to a file"
    outstring = list(zip(data['spacing'],\
                         data['error_2'],\
                         order_2,\
                         data['error_inf'],\
                         order_inf))
    header = ['#spacing', 'error_2 ', 'order_2 ', 'error_inf ', 'order_inf']
    # Format on the rows (: accepts the argument, < flushes left,
    # 20 denotes character width, .10e denotes scientific notation with
    # 10 in precision)
    header_format = "{:<20}" * (len(header))
    number_format = "{:<20.10e}" * (len(header))
    # * in front of header unpacks
    header_string = header_format.format(*header)
    text = header_string
    for string in outstring:
        text += '\n' + number_format.format(*string)
    print('\nNow printing the results of the convergence test:')
    print(text)
    # Write the found orders
    with open(os.path.join(root_folder, name + '.txt'), 'w' ) as f:
        f.write(text)
    print('\n')
#}}}

#{{{do_plot
def do_plot(data, order_2, order_inf, root_folder, name, extension,\
            xz_error_plot, xy_error_plot, show_plot, y_plane, z_plane,\
            directions):
    """Function which handles the actual plotting"""

    # Plot errors
    # Set the plotting style
    title_size = 30
    plt.rc("font", size = 30)
    plt.rc("axes", labelsize = 25, titlesize = title_size)
    plt.rc("xtick", labelsize = 25)
    plt.rc("ytick", labelsize = 25)
    plt.rc("legend", fontsize = 30)
    plt.rc("lines", linewidth = 3.0)
    plt.rc("lines", markersize = 20.0)
    plt_size = (10, 7)
    fig_no = 1
    # Try to make a figure with the current backend
    try:
        fig = plt.figure(fig_no, figsize = plt_size)
    except:
        # Switch if a backend needs the display
        plt.switch_backend('Agg')
        fig = plt.figure(fig_no, figsize = plt_size)

    ax = fig.add_subplot(111)

    # Plot errore
    # Plot the error-space plot for the 2-norm
    ax.plot(data['spacing'], data['error_2'], 'b-o', label=r'$L_2$')
    # Plot the error-space plot for the inf-norm
    ax.plot(data['spacing'], data['error_inf'], 'r-^', label=r'$L_\infty$')

    # Plot the order
    #{{{ Explanaition of the calculation
    # In the log-log plot, we have
    # ln(y) = a*ln(x) + ln(b)
    # y = error
    # x = spacing (found from linear regression)
    # a = order
    # b = where the straight line intersects with the ordinate
    #
    # Using the logarithmic identities
    # ln(x^a) = a*ln(x)
    # ln(x*y) = ln(x) + ln(y)
    # ln(x/y) = ln(x) - ln(y)
    #
    # We usually find b for x = 0. Here, on the other hand, we find it
    # by solving the equation for the smallest grid point:
    # ln[y(x[-1])] = a*ln(x[-1]) + ln(b),
    # so
    # ln(b) = ln[y(x[-1])] - a*ln(x[-1])
    # =>
    # ln(y) = a*ln(x) - a*ln(x[-1]) + ln[y(x[-1])]
    #       = a*[ln(x)-ln(x[-1])] + ln[y(x[-1])]
    #       = a*[ln(x/x[-1])] + ln[ys(x[-1])]
    #       = ln[(x/x[-1])^a*y(x[-1])]
    #}}}
    # Order in the inf norm
    ax.plot(\
             (data['spacing'][-1],\
              data['spacing'][0]),\
             (\
              ((data['spacing'][-1] / data['spacing'][-1])**order_inf)*\
                data['error_inf'][-1],\
              ((data['spacing'][0]  / data['spacing'][-1])**order_inf)*\
                data['error_inf'][-1]\
             ),\
             'm--',\
             label=r"$\mathcal{O}_{L_\infty}"+"%.2f"%(order_inf)+r"$")
    # Order in the 2 norm
    ax.plot(\
             (data['spacing'][-1],\
              data['spacing'][0]),\
             (\
              ((data['spacing'][-1] / data['spacing'][-1])**order_2)*\
                data['error_2'][-1],\
              ((data['spacing'][0]  / data['spacing'][-1])**order_2)*\
                data['error_2'][-1]\
             ),\
             'c--',\
             label=r"$\mathcal{O}_{L_2}"+"%.2f"%(order_2)+r"$")

    # Set logaraithmic scale
    ax.set_yscale('log')
    ax.set_xscale('log')

    # Set axis label
    ax.set_xlabel("Mesh spacing")
    ax.set_ylabel("Error norm")

    # Make the plot look nice
    # Plot the legend
    leg = ax.legend(loc="best", fancybox = True, numpoints=1)
    leg.get_frame().set_alpha(0.5)
    # Plot the grid
    ax.grid()
    # Make sure no collision between the ticks
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    # Includes the xlabel if outside
    plt.tight_layout()

    # Save the plot
    filename = os.path.join(root_folder, name  + '.' + extension)
    plt.savefig(filename)
    print('\nPlot saved to ' + filename + '\n'*2)

    if xz_error_plot:
        plot_xz_errors(data, root_folder, extension, y_plane, directions)

    if xy_error_plot:
        plot_xy_errors(data, root_folder, extension, z_plane, directions)

    if show_plot:
        plt.show()

    plt.close("all")
#}}}

#{{{plot_xz_errors
def plot_xz_errors(data, root_folder, extension, y_plane, directions):
    for E in data['error_field']:
        # Get mesh
        theta = E['dz'] * np.array(range(E['field'].shape[-1]))
        # We will plot without the ghost cells
        rho   = (E['dx'] * np.array(np.arange(0.5, E['field'].shape[-3])))
        THETA, RHO = np.meshgrid(theta, rho)

        # Check that we have only one ny
        if E['field'].shape[-2] != 1:
            if y_plane is None:
                y_plane = int(np.ceil(E['field'].shape[-2]/2) - 1)
                print('WARNING:')
                print('{} y planes found, plotting index {} (excluding ghost)'.\
                      format(E['field'].shape[-2], y_plane))
            E['field'] = E['field'][:, y_plane, :]

        if len(E['field'].shape) == 3:
            # Reshape the field (as we only have one y-plane)
            E['field'] =\
                    E['field'].reshape(E['field'].shape[-3], E['field'].shape[-1])

        fig = plt.figure(figsize = (10, 7))
        ax  = plt.subplot()
        cplot = ax.contourf(RHO, THETA, E['field'], 500, cmap=cm.RdYlBu_r)

        # Decorations
        cbar = plt.colorbar(cplot)
        cbar.ax.set_ylabel('Error')
        # Plot the grid
        ax.grid()
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel(r"$\theta$")
        # Set figure name
        name = ''
        if directions['dx']:
            name += 'nx={} '.format(E['field'].shape[0])
        if directions['dz']:
            name += 'nz={} '.format(E['field'].shape[1])
        fig.canvas.set_window_title(name)

        plt.yticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                           ['$0$', r'$\pi/2$', r'$\pi$',
                               r'$3\pi/2$', r'$2\pi$'])

        plt.tight_layout()

        # Save the plot
        filename = os.path.join(root_folder, name  + '.' + extension)
        plt.savefig(filename)
        print('\nPlot saved to ' + filename + '\n'*2)
#}}}

#{{{plot_xy_errors
def plot_xy_errors(data, root_folder, extension, z_plane, directions):
    for E in data['error_field']:
        # Get mesh
        z = E['dy'] * np.array(range(E['field'].shape[-2]))
        # We will plot without the ghost cells
        rho   = (E['dx'] * np.array(np.arange(0.5, E['field'].shape[-3])))
        Z, RHO = np.meshgrid(z, rho)

        # Check that we have only one nz
        if E['field'].shape[-1] != 1:
            if z_plane is None:
                z_plane = int(np.ceil(E['field'].shape[-1]/2) - 1)
                print('WARNING:')
                print('{} z planes found, plotting index {} (excluding ghost)'.\
                      format(E['field'].shape[-1], z_plane))
            E['field'] = E['field'][:, :, z_plane]

        if len(E['field'].shape) == 3:
            # Reshape the field (as we only have one y-plane)
            E['field'] =\
                    E['field'].reshape(E['field'].shape[-3], E['field'].shape[-2])

        fig = plt.figure(figsize = (10, 7))
        ax  = plt.subplot()
        cplot = ax.contourf(RHO, Z, E['field'], 500, cmap=cm.RdYlBu_r)

        # Decorations
        cbar = plt.colorbar(cplot)
        cbar.ax.set_ylabel('Error')
        # Plot the grid
        ax.grid()
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel(r"$z$")
        # Set figure name
        name = ''
        if directions['dx']:
            name += 'nx={} '.format(E['field'].shape[0])
        if directions['dy']:
            name += 'ny={} '.format(E['field'].shape[1])
        fig.canvas.set_window_title(name)

        plt.tight_layout()

        # Save the plot
        filename = os.path.join(root_folder, name  + '.' + extension)
        plt.savefig(filename)
        print('\nPlot saved to ' + filename + '\n'*2)
#}}}
