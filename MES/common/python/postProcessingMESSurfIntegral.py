#!/usr/bin/env python

"""Post processing which performs MES for the surface routines"""

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
                     extension     = 'png',\
                     show_plot     = False,\
                     use_dx        = False,\
                     use_dy        = False,\
                     use_dz        = False,\
                     ):
    """Collects the data members belonging to a convergence plot"""

    # Figure out the directions
    directions = {'dx':use_dx, 'dy':use_dy, 'dz':use_dz}

    # Make a variable to store the errors and the spacing
    data = {'e_Xout' :[],
            'e_Ydown':[],
            'e_Yup'  :[],
            'o_Xout' :[],
            'o_Ydown':[],
            'o_Yup'  :[],
            'spacing':[]}

    # Loop over the runs in order to collect
    for path in paths:
        # Collect the error
        e_Xout = collect('e_Xout' , path=path, info=False)
        e_Ydown= collect('e_Ydown', path=path, info=False)
        e_Yup  = collect('e_Yup'  , path=path, info=False)

        # Add the error field to data
        Lx = collect('Lx', path=path, info=False)
        Ly = collect('Ly', path=path, info=False)
        dx = collect('dx', path=path, info=False)
        dy = collect('dy', path=path, info=False)
        dz = collect('dz', path=path, info=False)

        # We are only interested in the absolute value
        data['e_Xout' ].append(np.abs(e_Xout ))
        data['e_Ydown'].append(np.abs(e_Ydown))
        data['e_Yup'  ].append(np.abs(e_Yup  ))

        # We are interested in the max of the spacings, so we make
        # an appendable list
        max_spacings = []
        # Collect the spacings
        if use_dx:
            max_spacings.append(np.max(dx))
        if use_dy:
            max_spacings.append(np.max(dy))
        if use_dz:
            max_spacings.append(np.max(dz))

        # Store the spacing in the data
        data['spacing'].append(np.max(max_spacings))

    # Sort the data
    data = sort_data(data)

    # Find the order of convergence in the 2 norm and infinity norm
    data = get_order(data)

    # Get the root name of the path (used for saving files)
    root_folder = os.path.join(paths[0].split('/')[0], "MES_results")
    # Set the top folder
    top_folder = []
    if use_dx:
        top_folder.append('dx')
    if use_dy:
        top_folder.append('dy')
    if use_dz:
        top_folder.append('dz')
    top_folder = "_".join(top_folder)
    # Combine
    full_folder = os.path.join(root_folder, top_folder)
    # Create folder if it doesn't exist
    if not os.path.exists(full_folder):
        os.makedirs(full_folder)
        print(full_folder + " created\n")

    # Get the name of the plot based on the first folder name
    name = paths[0].split('/')
    # Remove the root folder and put 'MES' in front
    name = 'MES'

    # Print the convergence rate
    print_convergence_rate(data, full_folder, name)

    # Plot
    do_plot(data, full_folder, name, extension, show_plot, directions)
#}}}

# Help functions
#{{{sort_data
def sort_data(data):
    """Sorts the data after highest grid spacing"""

    spacing_before_sort = data['spacing']

    for error in ['e_Xout', 'e_Ydown', 'e_Yup']:
        # Sort the data in case it is unsorted
        list_of_tuples_to_be_sorted =\
                list(zip(spacing_before_sort, data[error]))

        # Sort the list
        # Note that we are sorting in reverse order, as we want the
        # highest grid spacing first
        sorted_list = sorted(list_of_tuples_to_be_sorted, reverse = True)
        # Unzip the sorted list
        data['spacing'], data[error] = list(zip(*sorted_list))

    return data
#}}}

#{{{get_order
def get_order(data):

    for direction in ['Xout', 'Ydown', 'Yup']:
        # Initialize the orders
        data['o_'+direction].append(np.nan)

        # The order will be found by finding a linear fit between two
        # nearby points in the error-spacing plot. Hence, we must let
        # the index in the for loop run to the length minus one
        for index in range(len(data['spacing']) - 1):
            # p = polyfit(x,y,n) finds the coefficients of a polynomial p(x)
            # of degree that fits the data, p(x(i)) to y(i), in a least squares
            # sense.
            # The result p is a row vector of length n+1 containing the
            # polynomial coefficients in descending powers
            spacing_start   = np.log(data['spacing']     [index])
            spacing_end     = np.log(data['spacing']     [index + 1])
            error_start     = np.log(data['e_'+direction][index])
            error_end       = np.log(data['e_'+direction][index + 1])
            # Finding the order in the two norm
            polyfit_order = np.polyfit([spacing_start, spacing_end],\
                                       [error_start, error_end], 1)
            # Append it to the order
            data['o_'+direction].append(polyfit_order[0])

    return data
#}}}

#{{{print_convergence_rate
def print_convergence_rate(data, full_folder, name):
    "Prints the convergence rates to the screen and to a file"


    if "dx" in full_folder:
        dirs = ['Ydown', 'Yup']
    if "dy" in full_folder:
        dirs = ['Xout']
    if "dz" in full_folder:
        dirs = ['Xout', 'Ydown', 'Yup']

    for direction in dirs:
        outstring = list(zip(data['spacing'],\
                             data['e_'+direction],\
                             data['o_'+direction],\
                             ))
        header = ['#spacing', 'error', 'order']
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
        print('\nNow printing the results of the convergence test of {}:'\
                        .format(direction))
        print(text)
        # Write the found orders
        with open(os.path.join(full_folder,\
                  '{}_{}.txt'.format(name,direction)), 'w') as f:
            f.write(text)
        print('\n')
#}}}

#{{{do_plot
def do_plot(data, full_folder, name, extension, show_plot, directions):
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

    # Plot errors
    # Plot the error-space plot
    #{{{dx
    if "dx" in full_folder:
        ax.plot(data['spacing'], data['e_Ydown'], 'b-o',\
                label=r'Error lower plate')
        ax.plot(data['spacing'], data['e_Yup'], 'r-^',\
                label=r'Error upper plate')

        # Plot the order (see explanaition in postProcessingMES.py)
        ax.plot(\
                 (data['spacing'][-1],\
                  data['spacing'][0]),\
                 (\
                  ((data['spacing'][-1] / data['spacing'][-1])**data['o_Ydown'][-1])*\
                    data['e_Ydown'][-1],\
                  ((data['spacing'][0]  / data['spacing'][-1])**data['o_Ydown'][-1])*\
                    data['e_Ydown'][-1]\
                 ),\
                 'c--',\
                 label=r"$\mathcal{O}_{\mathrm{lower}}("+\
                        "%.2f"%(data['o_Ydown'][-1])+r")$")
        ax.plot(\
                 (data['spacing'][-1],\
                  data['spacing'][0]),\
                 (\
                  ((data['spacing'][-1] / data['spacing'][-1])**data['o_Yup'][-1])*\
                    data['e_Yup'  ][-1],\
                  ((data['spacing'][0]  / data['spacing'][-1])**data['o_Yup'][-1])*\
                    data['e_Yup'  ][-1]\
                 ),\
                 'm--',\
                 label=r"$\mathcal{O}_{\mathrm{upper}}("+\
                        "%.2f"%(data['o_Yup'][-1])+r")$")
    #}}}

    #{{{dy
    if "dy" in full_folder:
        ax.plot(data['spacing'], data['e_Xout'], 'b-o',\
                label=r'$\mathrm{Error outer plate}$')

        # Plot the order (see explanaition in postProcessingMES.py)
        ax.plot(\
                 (data['spacing'][-1],\
                  data['spacing'][0]),\
                 (\
                  ((data['spacing'][-1] / data['spacing'][-1])**data['o_Xout'][-1])*\
                    data['e_Xout'][-1],\
                  ((data['spacing'][0]  / data['spacing'][-1])**data['o_Xout'][-1])*\
                    data['e_Xout'][-1]\
                 ),\
                 'c--',\
                 label=r"$\mathcal{O}_{\mathrm{outer}}("+\
                        "%.2f"%(data['o_Xout'][-1])+r")$")
    #}}}

    #{{{dz
    if "dz" in full_folder:
        ax.plot(data['spacing'], data['e_Ydown'], 'b-o',\
                label=r'Error lower plate')
        ax.plot(data['spacing'], data['e_Yup'], 'r-^',\
                label=r'Error upper plate')
        ax.plot(data['spacing'], data['e_Xout'], 'k-v',\
                label=r'Error outer plate')

        # Plot the order (see explanaition in postProcessingMES.py)
        ax.plot(\
                 (data['spacing'][-1],\
                  data['spacing'][0]),\
                 (\
                  ((data['spacing'][-1] / data['spacing'][-1])**data['o_Ydown'][-1])*\
                    data['e_Ydown'][-1],\
                  ((data['spacing'][0]  / data['spacing'][-1])**data['o_Ydown'][-1])*\
                    data['e_Ydown'][-1]\
                 ),\
                 'c--',\
                 label=r"$\mathcal{O}_{\mathrm{lower}}("+\
                        "%.2f"%(data['o_Ydown'][-1])+r")$")
        ax.plot(\
                 (data['spacing'][-1],\
                  data['spacing'][0]),\
                 (\
                  ((data['spacing'][-1] / data['spacing'][-1])**data['o_Yup'][-1])*\
                    data['e_Yup'  ][-1],\
                  ((data['spacing'][0]  / data['spacing'][-1])**data['o_Yup'][-1])*\
                    data['e_Yup'  ][-1]\
                 ),\
                 'm--',\
                 label=r"$\mathcal{O}_{\mathrm{upper}}("+\
                        "%.2f"%(data['o_Yup'][-1])+r")$")

        ax.plot(\
                 (data['spacing'][-1],\
                  data['spacing'][0]),\
                 (\
                  ((data['spacing'][-1] / data['spacing'][-1])**data['o_Xout'][-1])*\
                    data['e_Xout'][-1],\
                  ((data['spacing'][0]  / data['spacing'][-1])**data['o_Xout'][-1])*\
                    data['e_Xout'][-1]\
                 ),\
                 'g--',\
                 label=r"$\mathcal{O}_{\mathrm{outer}}("+\
                        "%.2f"%(data['o_Xout'][-1])+r")$")
    #}}}
    # Set logaraithmic scale
    ax.set_yscale('log')
    ax.set_xscale('log')

    # Set axis label
    xlabel = []
    if directions['dx']:
        xlabel.append(r'$\Delta x$')
    if directions['dy']:
        xlabel.append(r'$\Delta y$')
    if directions['dz']:
        xlabel.append(r'$\Delta z$')
    xlabel = ", ".join(xlabel)
    xlabel = r'$\max ($' + xlabel + r'$)$'

    ax.set_xlabel(xlabel)
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
    filename = os.path.join(full_folder, name  + '.' + extension)
    plt.savefig(filename)
    print('\nPlot saved to ' + filename + '\n'*2)

    if show_plot:
        plt.show()

    plt.close("all")
#}}}
