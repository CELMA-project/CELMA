#!/usr/bin/env python

from boutdata import collect
import matplotlib.pyplot as plt
from pylab import plot
import numpy as np

from boututils.datafile import DataFile


"""Used to quick-check variables"""
def c2D(var_name=None, path='.', the_file = "BOUT.dmp.0.nc", prefix='BOUT.dmp',\
        xind=None, yind=None, zind=None, tind=None):
    print("Following variables are available:")
    with DataFile(path+'/'+the_file) as d:
        print(d.list())

    if prefix == 'BOUT.dmp':
        # Collect
        var = collect(var_name, path=path,\
                      xind=xind, yind=yind, zind=zind, tind=tind,\
                      yguards=True, xguards=True,\
                      info=False, prefix=prefix)
    else:
        # Read
        if xind == None:
            xind = [None, None]
        if yind == None:
            yind = [None, None]
        if zind == None:
            zind = [None, None]
        with DataFile(path+'/'+the_file) as d:
            dims = d.ndims(var_name)
            if dims == 0:
                ranges = None
            elif dims == 1:
                ranges = [xind[0], xind[1]]
            elif dims == 2:
                ranges = [xind[0], xind[1], yind[0], yind[1]]
            elif dims == 3:
                ranges = [xind[0], xind[1], yind[0], yind[1], zind[0], zind[1]]
            var = d.read(var_name, ranges=ranges)

    if hasattr(var, "__iter__"):
        if len(var.shape) == 4:
            if var.shape[0] == 1:
                # Recast to three dimensions
                new_var = np.empty(var.shape[1:None])
                new_var = var[0,:,:,:]
                var = new_var
            # Check if we are operating with 2D variables
            if var.shape[0] == 1:
                # Recast to 2 dimensions
                new_var = np.empty(var.shape[1:None])
                new_var = var[0,:,:]
                var = new_var
                # Check if 1D
                if var.shape[0] == 1:
                    new_var = np.empty(var.shape[1:None])
                    new_var = var[0,:]
                    var = new_var
                    # Check if 0D
                    if var.shape[0] == 1:
                        var = var[0]
            elif var.shape[1] == 1:
                # Recast to 2 dimensions
                new_var = np.empty(var.shape[1:None])
                new_var = var[:,0,:]
                var = new_var
                # Check if 1D
                if var.shape[1] == 1:
                    new_var = np.empty(var.shape[1:None])
                    new_var = var[:,1]
                    var = new_var
            elif var.shape[2] == 1:
                # Recast to 2 dimensions
                new_var = np.empty(var.shape[1:None])
                new_var = var[:,:,0]
                var = new_var

        if len(var.shape) == 2:
            # Start fig
            fig = plt.figure()
            ax = fig.add_subplot(111)
            the_plot = ax.imshow(var.transpose(), interpolation='none')
            ax.set_xlim = (0, var.shape[0]-1)
            ax.set_ylim = (0, var.shape[1]-1)

            # Add colorbar
            plt.colorbar(the_plot)
            # Plot
            plt.show()
    else:
        print(var_name + "=" + str(var))

    return var
