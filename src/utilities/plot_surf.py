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
    :return data: np.array of values
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

def get_beam_position(X,Y,plane):
    '''
    Routine to find bounding points of beam
    :param X:
    :param Y:
    :return:
    '''
    xmin = 10
    ymin = 10
    if plane == 'y':
        for i in range(len(X)):
            if (X[i] < 0) and (abs(X[i]) < abs(xmin)) and (abs(Y[i]) < 0.01 ):
                xmin = X[i]
        for i in range(len(Y)):
            if (Y[i] < 0) and (abs(Y[i]) < abs(ymin)) and (abs(X[i]) < 0.05):
                ymin = Y[i]
        width = 1
        height = 1

    elif plane == 'z':
        for i in np.nditer(X):
            if (i < 0) and (abs(i) < abs(xmin)):
                xmin = i
        width = 1
        height = 4

        ymin = 0

    dict = {'xmin':xmin,'ymin':ymin,'width':width,'height':height}
    return dict

def plot_instantaneous(X,Y,Z,figname,figuresize=(7.54,2.5),basic=False,colormap='rainbow',plane='y'):
    '''
    Routine plots the instantaneous values for a parameter
    of choosing
    :param X:
    :param Y:
    :param Z:
    :param figname: Name of output file
    :return:
    '''
    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    # Create figure and axes
    fig,ax = plt.subplots(figsize=figuresize)

    # Plot the points
    #plt.scatter(x_c, y_c, marker='o', s=5, zorder=10)

    ##############################
    # BASIC PLOT or FULL PLOT
    ##############################

    # Plot the basic data, similar to cell plots in paraview
    if(basic):
        zic = griddata((X, Y), Z, (X[None, :], Y[:, None]))
        plt.pcolormesh(X, Y, zic, cmap=plt.get_cmap(colormap))
    else:
        # Plot interpolated data
        # Create grid for interpolations
        # Last input parameter is the number of grid points (10000) is large
        xic = np.linspace(min(X), max(X), 1000)
        yic = np.linspace(min(Y), max(Y), 1000)

        ##############################
        bands = 20
        ##############################

        zic = griddata((X, Y), Z, (xic[None,:],yic[:,None]))

        plt.contourf(xic,yic,zic,bands,cmap=plt.get_cmap(colormap))

        cbar = plt.colorbar()
        #cbar.set_clim(vmin=-1, vmax=0.4)
        #cbar.set_ticks([-1,-0.5,0,0.5])
        cbar.set_label("Pressure")

        C = plt.contour(xic, yic, zic, 20, colors='black')
        #plt.clabel(C, inline=1, fontsize=10)



        #plt.contour(xic,yic,zic,12,linewidths = 0.5, colors = 'k')
        #plt.pcolormesh(xic, yic, zic,cmap = plt.get_cmap('spectral'))

    # Create a Rectangle patch
    # Find position of beam
    minimums = get_beam_position(X,Y,plane)
    xmin = minimums['xmin'];ymin = minimums['ymin']
    wdith = minimums['width'];height = minimums['height']
    plt.ylim((-4,4))
    #plt.xlim((-3, 15))

    # Add the patch to the Axes
    rect = patches.Rectangle((xmin,ymin),wdith,height,edgecolor='grey',facecolor='grey',alpha=1.0,zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.close()
    #plt.show()

def plot_average(X,Y,Z,figname,figuresize=(7.54,2.5),basic=False,colormap='rainbow',plane='y'):
    '''
        Routine plots the instantaneous values for a parameter
        of choosing
        :param X:
        :param Y:
        :param Z:
        :param figname: Name of output file
        :return:
        '''
    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    # Create figure and axes
    fig, ax = plt.subplots(figsize=figuresize)
    plt.contourf(X, Y, Z, 20, cmap=plt.get_cmap(colormap))

    cbar = plt.colorbar()
    # cbar.set_clim(vmin=-1, vmax=0.4)
    # cbar.set_ticks([-1,-0.5,0,0.5])
    cbar.set_label("Pressure")

    C = plt.contour(X, Y, Z, 20, colors='black')
    # plt.clabel(C, inline=1, fontsize=10)

    # Create a Rectangle patch
    # Find position of beam
    plt.ylim((-4, 4))
    # plt.xlim((-3, 15))

    # Add the patch to the Axes
    rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=1.0, zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname, bbox_inches='tight', dpi=300)
    plt.close()
    # plt.show()

def get_data(infile,parallel,value,plane,npes):
    '''
    Routine to get the data from a file name
    :param fname:
    :param parallel:
    :return:
    '''
    if plane == 'y':
        Y=2
    elif plane == 'z':
        Y=1
    # Extract Data for given timestep
    if (parallel):
        x_c = []
        y_c = []
        z_c = []
        for i in range(npes):
            fname = 'processor' + str(i) + '/postProcessing/' + infile
            if os.path.isfile(fname):
                data = read_csv(fname)
                x_c.extend(data[:, 0])
                y_c.extend(data[:, Y])
                z_c.extend(data[:, value])
    else:
        fname = 'postProcessing/' + infile
        if os.path.isfile(fname):
            data = read_csv(fname)
            x_c = data[:, 0]
            y_c = data[:, Y]
            z_c = data[:, value]

    x_c = np.ravel(np.array(x_c))
    y_c = np.ravel(np.array(y_c))
    z_c = np.ravel(np.array(z_c))

    dict = {'X':x_c,'Y':y_c,'Z':z_c}
    return dict

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

    ##############################

    csvFile = "yplane_4_119.625.csv"
    outFile = "yplane_4_119.pdf"

    # 3 - Pressure
    # (4,5,6) - Ux,Uy,Uz
    property = 3
    ##############################

    print("\n-------------------------------------------")
    print("Plotting Instantaneous Value, file: {0}".format(csvFile))
    print("-------------------------------------------\n")

    # Get the data
    data = get_data(csvFile,parallel,3,'y',npes)

    # Plot the Instantaneous flow
    plot_instantaneous(data['X'],data['Y'],data['Z'],outFile,plane='y')

    #####################################################################################

    ##############################
    plane_name ="yplane_4"
    csvFile = "yplane_4_118.5.csv"
    outFile = "yplane_4_average.pdf"
    ##############################

    print("\n-------------------------------------------")
    print("Plotting Average for {0}".format(plane_name))
    print("-------------------------------------------\n")


    # Find number of timeSteps in a dir and get a sorted file List
    for i in range(npes):
        fname = 'processor' + str(i) + '/postProcessing/' + csvFile
        if os.path.isfile(fname):
            path = 'processor' + str(i) + '/postProcessing/'
            # Create list of file numbers
            tmp_list = []
            for file in os.listdir(path):
                if plane_name in file:
                    tmp = file.replace(plane_name,'').replace('_','').replace('.csv','')
                    if tmp == '0':
                        continue
                    else:
                        tmp_list.append(tmp)
            exit
    file_list = sorted(tmp_list)
    time_count = len(file_list)

    print("Number of time steps: {0}".format(time_count))

    # CREATE BASE GRID
    # Create base grid with which to map the new values onto
    # Load in a single time
    first = True
    for i in file_list:
        fname=plane_name+'_'+str(i)+'.csv'
        data = get_data(fname, parallel, 3, 'y', npes)
        if first:
            # CREATE BASE GRID
            xic = np.linspace(min(data['X']), max(data['X']), 1000)
            yic = np.linspace(min(data['Y']), max(data['Y']), 1000)
            first = False
            zic_average = griddata((data['X'], data['Y']), data['Z'], (xic[None, :], yic[:, None]))
        else:
            zic_average = np.add(zic_average,griddata((data['X'], data['Y']),
                                data['Z'], (xic[None, :], yic[:, None])))

    zic_average = zic_average/time_count

    # Plot the Average flow
    plot_average(xic, yic, zic_average, outFile, plane='y')





