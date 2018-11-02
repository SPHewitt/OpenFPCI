#!/usr/bin/python2.7

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc

import numpy as np

from scipy.interpolate import griddata


def read_csv(fname):
    '''
    Routine reads a csv file with ',' delimiter
    :param fname: The filename to be read
    :return: data: np.array of values
    '''
    data = np.loadtxt(fname,delimiter=',')
    return data

def check_parallel():
    '''
    Very crude method of checking if case was run in serial
    or Parallel. It looks at the dir list and searches for a
    dir with the prefix processor
    :return: flag: bool
    '''
    flag=False

    dir_list = os.listdir(os.getcwd())
    for i in dir_list:
        if 'processor' in i:
            flag = True
            return flag

    return flag

def get_num_proc():
    '''
    Crude way of counting processors
    :return: npes: number of processors
    '''
    npes = 0

    dir_list = os.listdir(os.getcwd())
    for i in dir_list:
        if 'processor' in i:
            npes=npes+1

    return npes

if __name__ == '__main__':

    # Heading
    print("\n-------------------------------------------")
    print("Plotting Surface ")
    print("-------------------------------------------\n")

    # Check if parallel
    parallel = check_parallel()

    # Get Number of processors
    if (parallel):
        npes = get_num_proc()
        print("Number of Processors: {0}".format(npes))

    ###################
    csvFile = "zplane_0.1.csv"
    ##################

    # Extract Data for given timestep
    if(parallel):
        x_c = []
        y_c = []
        z_c = []
        for i in range(npes):
            fname='processor'+str(i)+'/postProcessing/'+csvFile
            if os.path.isfile(fname):
                data = read_csv(fname)
                x_c.extend(data[:, 0])
                y_c.extend(data[:, 1])
                z_c.extend(data[:, 4])
    else:
        fname='postProcessing/'+csvFile
        if os.path.isfile(fname):
            data = read_csv(fname)
            x_c = data[:,0]
            y_c = data[:,1]
            z_c = data[:,4]

    x_c=np.ravel(np.array(x_c))
    y_c=np.ravel(np.array(y_c))
    z_c=np.ravel(np.array(z_c))

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    # Create figure and axes
    fig,ax = plt.subplots(figsize=(5.51,3.5))

    # Plot the points
    #plt.scatter(x_c, y_c, marker='o', s=5, zorder=10)

    # Plot the basic data, similar to cell plots in paraview
    #zic = griddata((x_c, y_c), z_c, (x_c[None, :], y_c[:, None]))
    #plt.pcolormesh(x_c, y_c, zic, cmap=plt.get_cmap('rainbow'))

    # Plot interpolated data
    # Create grid for interpolations
    # Last input parameter is the number of grid points (10000) is large
    xic = np.linspace(min(x_c), max(x_c), 1000)
    yic = np.linspace(min(y_c), max(y_c), 1000)

    zic = griddata((x_c, y_c), z_c, (xic[None,:],yic[:,None]))
    plt.contourf(xic,yic,zic,12,cmap=plt.get_cmap('RdGy'),linewidths = 0.5)
    cbar = plt.colorbar()
    cbar.set_label("Pressure")
    plt.contour(xic, yic, zic, 12,linewidths=0.5,colors='k')

    #plt.contour(xic,yic,zic,12,linewidths = 0.5, colors = 'k')
    #plt.pcolormesh(xic, yic, zic,cmap = plt.get_cmap('spectral'))



    # Create a Rectangle patch
    # Add the patch to the Axes
    #rect = patches.Rectangle((-0.5,-0.5),1,1,linewidth=1,edgecolor='k',facecolor='k',alpha=1.0)
    rect = patches.Rectangle((-0.5, 0), 1, 4, linewidth=1, edgecolor='k', facecolor='k', alpha=1.0)
    ax.add_patch(rect)


    plt.show()
