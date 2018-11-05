#!/usr/bin/python2.7

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc

import numpy as np

from scipy.interpolate import griddata

#### Misc Functions
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

def create_average(dict,plane_name):
    '''
    Routine creates a base grid based and interpolates all
    files in file list onto gird, and averages
    :param dict:
    :param plane_name:
    :return:
    '''
    # CREATE BASE GRID
    # Create base grid with which to map the new values onto
    # Load in a single time
    file_list=dict['files']
    time_count=dict['time_count']

    first = True
    for i in file_list:
        fname = plane_name + '_' + str(i) + '.csv'
        data = get_data(fname, parallel, npes)
        if first:
            # CREATE BASE GRID
            xic = np.linspace(min(data['X']), max(data['X']), 1000)
            yic = np.linspace(min(data['Z']), max(data['Z']), 1000)
            first = False
            zic_average = griddata((data['X'], data['Z']), data[value], (xic[None, :], yic[:, None]))
        else:
            zic_average = np.add(zic_average, griddata((data['X'], data['Z']),
                                                       data[value], (xic[None, :], yic[:, None])))

    zic_average = zic_average / time_count
    out_dict={'xic':xic,'yic':yic,'zic_average':zic_average}
    return out_dict

#### Get Functions ####
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

def get_data(infile,parallel,npes):
    '''
    Routine to get the data from a file name
    :param fname:
    :param parallel:
    :return:
    '''
    # Extract Data for given timestep
    if (parallel):
        x_c = []
        y_c = []
        z_c = []
        p_c = []
        Ux_c = []
        Uy_c = []
        Uz_c = []
        for i in range(npes):
            fname = 'processor' + str(i) + '/postProcessing/' + infile
            if os.path.isfile(fname):
                data = read_csv(fname)
                x_c.extend(data[:, 0])
                y_c.extend(data[:, 1])
                z_c.extend(data[:, 2])
                p_c.extend(data[:, 3])
                Ux_c.extend(data[:, 4])
                Uy_c.extend(data[:, 5])
                Uz_c.extend(data[:, 6])
    else:
        fname = 'postProcessing/' + infile
        if os.path.isfile(fname):
            data = read_csv(fname)
            x_c = data[:, 0]
            y_c = data[:, 1]
            z_c = data[:, 2]
            p_c = data[:, 3]
            Ux_c = data[:, 4]
            Uy_c = data[:, 5]
            Uz_c = data[:, 6]

    x_c = np.ravel(np.array(x_c))
    y_c = np.ravel(np.array(y_c))
    z_c = np.ravel(np.array(z_c))
    p_c = np.ravel(np.array(p_c))
    Ux_c = np.ravel(np.array(Ux_c))
    Uy_c = np.ravel(np.array(Uy_c))
    Uz_c = np.ravel(np.array(Uz_c))

    dict = {'X':x_c,'Y':y_c,'Z':z_c,'pressure':p_c,'Ux':Ux_c,'Uy':Uy_c,'Uz':Uz_c}
    return dict

def get_file_list(plane_name,npes):
    '''
    Get a list of timestep files
    :return:
    '''
    # Find number of time steps in a dir and get a sorted file List
    flag = False
    for i in range(npes):
        path = 'processor' + str(i) + '/postProcessing/'
        for j in os.listdir(path):
            if plane_name in j:
                flag = True
        if (flag):
            break

    # Create list of file numbers
    tmp_list = []
    for file in os.listdir(path):
        if plane_name in file:
            tmp = file.replace(plane_name, '').replace('_', '').replace('.csv', '')
            if tmp == '0':
                continue
            else:
                tmp_list.append(tmp)
        exit

    file_list = sorted(tmp_list)
    time_count = len(file_list)
    dict={'files':file_list,'time_count':time_count}
    return dict

#### Plotting Functions ####
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

def plot_average(X,Y,Z,figname,figuresize=(7.54,2.5),basic=False,colormap='rainbow',plane='y',label='Colorbar'):
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
    cbar.set_label(label)

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
    ##############################

    print("\n-------------------------------------------")
    print("Plotting Instantaneous Value, file: {0}".format(csvFile))
    print("-------------------------------------------\n")

    # Get the data
    data = get_data(csvFile,parallel,npes)

    # Plot the Instantaneous flow
    plot_instantaneous(data['X'],data['Z'],data['pressure'],outFile,plane='y')

    #####################################################################################
    #####################################################################################
    #####################################################################################
    # PLOTTING AVERAGES

    ##############################
    plane_name ="yplane_4"
    value='Uy' #'pressure,Ux,Uy,Uz

    outFile = plane_name+"_average_"+value+".pdf"
    ##############################

    print("\n-------------------------------------------")
    print("Plotting Average for {0}".format(plane_name))
    print("-------------------------------------------\n")

    dict = get_file_list(plane_name, npes)

    print("Number of time steps: {0}".format(dict['time_count']))
    dict = create_average(dict, plane_name)
    xic = dict['xic']
    yic = dict['yic']
    zic_average = dict['zic_average']

    # Plot the Average flow
    plot_average(xic, yic, zic_average, outFile, plane='y',label=value)

    #####################################################################################
    #####################################################################################
    #####################################################################################

    sys.exit()
    # Plotting Vorticity
    data = get_data("yplane_4_119.625.csv", parallel, npes)

    dUx=np.gradient(data['Ux'])
    dUy=np.gradient(data['Uy'])
    dUz=np.gradient(data['Uz'])

    dx=np.gradient(data['X'])
    dy=np.gradient(data['Y'])
    dz=np.gradient(data['Z'])

    print data['Z']
    print(dz)

    #vorticity = (dUz/dy - dUy/dz) + (dUx/dz - dUz/dx) + (dUy/dx-dUx/dy)
    vorticity = np.divide(dUy,dz)



