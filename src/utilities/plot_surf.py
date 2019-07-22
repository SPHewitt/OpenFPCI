#!/usr/bin/python3.5

import os
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc
import numpy as np
from scipy.interpolate import griddata
import matplotlib as mpl


from scipy.interpolate import spline
import scipy.ndimage

# Property dictionary
property_dictionary={'x':0,'y':1,'z':2,'p':3,'pmean':4,'ux':5,'uy':6,'uz':7,'Ux':8,'Uy':9,'Uz':10,
                     'UU':11,'VV':12,'WW':13,'UV':14,'VW':15,'UW':16,'wx':17,'wy':18,'wz':19}


# Plane lookups
plane_lookup={'x':('Z','Y'),'y':('X','Z'),'z':('X','Y')}

# Global variable
PARALLEL = True


BINARY=True
# If Binary files
# if len(sys.argv) > 1:
#    BINARY = True
# else:
#    BINARY = False

print("Binary: {}".format(BINARY))

# Misc Functions
def read_csv(fname):

    '''
    Routine reads a csv file with ',' delimiter
    :param fname: The filename to be read
    :return data: np.array of values
    '''
    global BINARY

    if BINARY:
        data = np.fromfile(fname, dtype=np.float64)
        length = len(data)
        cols = 20
        rows = int(length/cols)
        data = data.reshape((rows,cols))
    else:
        data = np.loadtxt(fname,delimiter=',')
    return data


def read_bin(fname):
    # Declare a datatype
    my_dtype = np.dtype(
        [("x", np.float64),  ("y", np.float64), ("z", np.float64), ("p", np.float64), ("pMean", np.float64),
         ("ux", np.float64), ("uy", np.float64), ("uz", np.float64), ("Ux", np.float64), ("Uy", np.float64),
         ("Uz", np.float64), ("UU", np.float64), ("VV", np.float64), ("WW", np.float64), ("UV", np.float64),
         ("VW", np.float64), ("UW", np.float64), ("wx", np.float64), ("wy", np.float64), ("wz", np.float64)])

    # Read file
    data = np.fromfile(fname, dtype=my_dtype)

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


def create_average(dict,plane_name,prop):

    '''
    Routine creates a base grid based and interpolates all
    files in file list onto gird, and averages
    :param dict: Dictionary of files, keys = str(time) and values = float(time)
                e.g. {'90':90,'90.5',90.5}
    :param plane_name: name of the plane
    :return:
    '''

    global BINARY

    # CREATE BASE GRID
    # Create base grid with which to map the new values onto
    time_count=len(dict.keys())
    file_list=sorted(dict.keys(),key=float)

    # Plane lookup
    X_name,Y_name = plane_lookup[plane_name[0]]
    if plane_name[0] == 'x':
        Y_name, X_name = plane_lookup[plane_name[0]]

    first = True
    for i in file_list:
        if BINARY:
            fname = plane_name + '_' + str(i)
        else:
            fname = plane_name + '_' + str(i) + '.csv'
        data = get_data(fname, npes,prop=prop)
        if first:
            # CREATE BASE GRID
            xic = np.linspace(min(data[X_name]), max(data[X_name]), 1000)
            yic = np.linspace(min(data[Y_name]), max(data[Y_name]), 1000)
            zic = griddata((data[X_name], data[Y_name]), data['C'], (xic[None, :], yic[:, None]))
            first = False
        else:
            try:
                zic = np.add(zic, griddata((data[X_name], data[Y_name]),data['C'], (xic[None, :], yic[:, None])))
            except Exception:
                pass


    zic = np.divide(zic,time_count)
    out_dict={'xic':xic,'yic':yic,'zic':zic}

    return out_dict


def write_to_file(data,filename):
    # Write the array to disk
    with open(filename, 'w') as outfile:
        # I'm writing a header here just for the sake of readability
        # Any line starting with "#" will be ignored by numpy.loadtxt
        outfile.write('# Array shape: {0}\n'.format(data.shape))

        # Iterating through a ndimensional array produces slices along
        # the last axis. This is equivalent to data[i,:,:] in this case
        for data_slice in data:

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.
            np.savetxt(outfile, data_slice, fmt='%-7.2f')


def load_from_file(filename):
    # Read the array from disk
    new_data = np.loadtxt(filename)
    # However, going back to 3D is easy if we know the
    # original shape of the array
    new_data = new_data.reshape(1000,1000)
    return new_data


def get_limits(plane_name,value):
    if value == 'pmean':
        if plane_name[0] == 'x':
            xlim = (-4, 4)
            ylim = (0, 6)
        if plane_name[0] == 'y':
            xlim = (-4, 10)
            ylim = (-4, 4)
        if plane_name[0] == 'z':
            xlim = (-4, 10)
            ylim = (0, 6)

    elif value == 'wx':
        if plane_name[0] == 'x':
            xlim = (-3, 3)
            ylim = (0, 5)
        if plane_name[0] == 'y':
            xlim = (-4, 6)
            ylim = (-3, 3)
        if plane_name[0] == 'z':
            xlim = (-4, 10)
            ylim = (0, 6)

    else:
        if plane_name[0] == 'x':
            xlim = (-4, 4)
            ylim = (0, 6)
        if plane_name[0] == 'y':
            xlim = (-4, 10)
            ylim = (-4, 4)
        if plane_name[0] == 'z':
            xlim = (-4, 10)
            ylim = (0, 6)

    return xlim,ylim

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
            if (X[i] < 0) and (abs(X[i]) < abs(xmin)) and (abs(Y[i]) < 0.05 ):
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

    #print xmin
    #print ymin
    dict = {'xmin':xmin,'ymin':ymin,'width':width,'height':height}
    return dict


def get_data(infile,npes,prop):
    '''
    Routine to get the data from a file name
    :param fname:
    :param npes:
    :param prop:
    :return:
    Data formatted as follows:
    x y z p pmean ux uy uz Ux Uy Uz UU VV WW UV VW UW wx wy wz
    '''
    global PARALLEL

    # Get the value in the array
    value = property_dictionary[prop]


    # Extract Data for given timestep
    if (PARALLEL):
        x_c = []
        y_c = []
        z_c = []
        prop_c = []

        for i in range(npes):
            fname = 'processor' + str(i) + '/postProcessing/' + infile
            if os.path.isfile(fname):
                data = read_csv(fname)
                if data.ndim == 0:
                    continue
                elif data.ndim == 1:
                    x_c.append(data[0])
                    y_c.append(data[1])
                    z_c.append(data[2])
                    prop_c.append(data[value])
                elif data.ndim == 2:
                    x_c.extend(data[:, 0])
                    y_c.extend(data[:, 1])
                    z_c.extend(data[:, 2])
                    prop_c.extend(data[:, value])

    else:
        fname = 'postProcessing/' + infile
        if os.path.isfile(fname):
            data = read_csv(fname)
            x_c = data[:, 0]
            y_c = data[:, 1]
            z_c = data[:, 2]
            prop_c = data[:, value]

    x_c = np.ravel(np.array(x_c))
    y_c = np.ravel(np.array(y_c))
    z_c = np.ravel(np.array(z_c))
    prop_c = np.ravel(np.array(prop_c))

    # Error checking for empty arrays
    if len(x_c) == 0:
        print('No data in arrays')
        sys.exit()

    dict = {'X':x_c,'Y':y_c,'Z':z_c,'C':prop_c}
    return dict


def get_turbulent_data(infile,npes):
    '''
    Routine to get the data from a file name
    :param fname:
    :param npes:
    :param prop:
    :return:
    Data formatted as follows:
    x y z p pmean ux uy uz Ux Uy Uz UU VV WW UV VW UW wx wy wz
    '''
    global PARALLEL

    # Get the value in the array
    valueA = property_dictionary['UU']
    valueB = property_dictionary['Ux']

    # Extract Data for given timestep
    if (PARALLEL):
        x_c = []
        y_c = []
        z_c = []
        prop_c = []

        for i in range(npes):
            fname = 'processor' + str(i) + '/postProcessing/' + infile
            if os.path.isfile(fname):
                data = read_csv(fname)
                if data.ndim == 0:
                    continue
                elif data.ndim == 1:
                    x_c.append(data[0])
                    y_c.append(data[1])
                    z_c.append(data[2])
                    prop_c.append(data[valueA]-(np.square(data[valueB])))
                elif data.ndim == 2:
                    x_c.extend(data[:, 0])
                    y_c.extend(data[:, 1])
                    z_c.extend(data[:, 2])
                    prop_c.extend(data[:, valueA]-(np.square(data[:,valueB])))

    else:
        fname = 'postProcessing/' + infile
        if os.path.isfile(fname):
            data = read_csv(fname)
            x_c = data[:, 0]
            y_c = data[:, 1]
            z_c = data[:, 2]
            prop_c = data[:, value]

    x_c = np.ravel(np.array(x_c))
    y_c = np.ravel(np.array(y_c))
    z_c = np.ravel(np.array(z_c))
    prop_c = np.ravel(np.array(prop_c))

    # Error checking for empty arrays
    if len(x_c) == 0:
        print('No data in arrays')
        sys.exit()

    dict = {'X':x_c,'Y':y_c,'Z':z_c,'C':prop_c}
    return dict


def get_file_list(plane_name,npes,verbose=True):
    '''
    Find number of time steps in a dir and get a sorted file List
    create a dictionary with keys = to name of file and values = floats
    {'100':100}
    :return:
    '''

    dict = {}
    for i in range(npes):
        path = 'processor' + str(i) + '/postProcessing/'
        for j in os.listdir(path):
            if plane_name in j:
                tmp=j.replace(plane_name, '').replace('_', '').replace('.csv', '')
                if tmp == '0':
                    continue
                else:
                    dict.setdefault(tmp,float(tmp))

    if(verbose):
        print("Start Time: {0}".format(min(dict.values())))
        print("End Time: {0}".format(max(dict.values())))
        print("Total Steps: {0}".format(len(dict.keys())))

    return dict

#### Plotting Functions ####


def plot_basic(X,Y,Z,figname,figuresize=(7.54,2.5),colormap='rainbow'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    # Create figure and axes
    fig, ax = plt.subplots(figsize=figuresize)

    # Plot the points
    # plt.scatter(x_c, y_c, marker='o', s=5, zorder=10)

    zic = griddata((X, Y), Z, (X[None, :], Y[:, None]))
    plt.pcolormesh(X, Y, zic, cmap=plt.get_cmap(colormap))

    plt.savefig(figname,bbox_inches='tight',dpi=300)


def plot_instantaneous(data,figname,bands=100,c_label='cbar',figuresize=(7.54,2.5),colormap='jet',plane='y',v_min=-1,v_max=1,step=0.2,xlim=None,ylim=None):
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
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
    rc('text', usetex=True)

    # Plane lookup
    X_name,Y_name = plane_lookup[plane]
    if plane=='x':
        X = data[Y_name];Y = data[X_name];Z = data['C']
    else:
        X = data[X_name];Y = data[Y_name]; Z = data['C']


    # Create figure and axes
    fig,ax = plt.subplots(figsize=figuresize)

    # Create grid for interpolations
    xic = np.linspace(min(X), max(X), 1000)
    yic = np.linspace(min(Y), max(Y), 1000)

    zic = griddata((X, Y), Z, (xic[None,:],yic[:,None]))

    if plane == 'x':
        zic=zic.transpose()
        plt.contourf(yic,xic,zic,bands,cmap=plt.get_cmap(colormap))
        plt.contourf(yic, xic, zic, bands, cmap=plt.get_cmap(colormap))
    else:
        #plt.contourf(xic,yic,zic,bands,cmap=plt.get_cmap(colormap),vmin=v_min,vmax=v_max)
        #plt.contourf(xic, yic, zic, bands, cmap=plt.get_cmap(colormap), vmin=v_min, vmax=v_max)
        plt.contourf(xic,yic,zic,bands,cmap=plt.get_cmap(colormap))
        plt.contourf(xic, yic, zic, bands, cmap=plt.get_cmap(colormap))

    if xlim is None:
        pass
    else:
        plt.xlim=(xlim)

    if ylim is None:
        pass
    else:
        plt.ylim=(ylim)


    m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    m.set_clim(v_min, v_max)
    m.set_array(Z)
    cbarlabels = np.arange(v_min, v_max + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)
    cbar = plt.colorbar(m)
    cbar.set_ticks(cbarlabels)
    cbar.set_ticklabels(cbarlabels)
    cbar.set_label(c_label)
    cityC = plt.contour(xic, yic, zic, bands, linewidths=0.2, colors='black', alpha=0.2)
    # plt.clabel(C, inline=1, fontsize=10)

    # Create a Rectangle patch
    if plane == 'y':
        minimums = get_beam_position(X,Y,plane)
        xmin = minimums['xmin'];ymin = minimums['ymin']
        width = minimums['width'];height = minimums['height']
        rect = patches.Rectangle((xmin, ymin), width, height, edgecolor='grey', facecolor='grey', alpha=1.0, zorder=10)
        ax.add_patch(rect)
        plt.xlim((-3,5))
        plt.ylim((-2,2))

    elif plane == 'z':
        xmin = -0.5; ymin=0; width =1; height = 4
        rect = patches.Rectangle((xmin, ymin), width, height, edgecolor='grey', facecolor='grey',fill=False ,alpha=1.0, zorder=10)
        ax.add_patch(rect)

    else:
        rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', alpha=1.0, linestyle='--',fill=False)
        ax.add_patch(rect)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.close()


def plot_streamlines_zplane_100(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='rainbow',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',figname='zplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-5,10))
    plt.ylim((0,6))

    colors = plt.cm.jet(np.linspace(0,1,4))
    #plt.contourf(X, Y, P, 256,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    # Upstream
    x_values=np.arange(-5,-4,0.5)
    y_values=np.arange(0.1,6,0.5)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=15,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # Horseshoe
    x_values=np.arange(-2,-0.5,0.5)
    y_values=np.arange(0.005,0.5,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=10,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'go', markersize=1)

    # Near Wake
    x_values=np.arange(1,2,0.5)
    y_values=np.arange(0.1,4,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=10,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)

    # Far wake
    x_values=np.arange(7,8,1)
    y_values=np.arange(0.1,4,0.2)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=5,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ko',markersize=1)

    # Added to obtain saddle point for re=100
    x_values=np.arange(4,6,0.25)
    y_values=np.arange(1,3,0.25)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=10,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)

    # Added to obtain saddle point for re=100
    x_values=np.full(1,5.2)
    y_values=np.arange(1.3,1.75,0.2)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)

    # Added to obtain saddle point for re=100
    x_values=np.full(1,5.37)
    y_values=np.arange(1.47,1.75,0.12)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)


    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=2,maxlength=100.0,color='black',arrowsize=0.4)

    plt.plot(1.64,3.4,color='purple',marker='P',markersize=8)
    plt.text(1.8,3.6,'D',color='purple',fontsize=16,weight='black')
    #plt.plot(0.714,2.29, color='purple',marker='P',markersize=8)
    plt.plot(5.31, 1.59, color='purple',marker='P',markersize=8)
    plt.text(5.5,1.8,'I',color='purple',fontsize=16,weight='black')
    #plt.plot(-1.06, 0.203, color='purple',marker='P',markersize=5)

    #plt.text(1.64,3.4,'D',color='purple')
    #plt.text(5.31, 1.59,'I' ,color='purple')

    rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('zplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)

    plt.close()


def plot_streamlines_zplane_250(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='rainbow',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',figname='zplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-5,10))
    plt.ylim((0,6))

    colors = plt.cm.jet(np.linspace(0,1,4))
    #plt.contourf(X, Y, P, 256,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    myColor = 'dimgrey'

    # Upstream
    x_values=np.arange(-5,-4,0.5)
    y_values=np.arange(0.1,6,0.5)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color=myColor,arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'yo', markersize=1)

    # Horseshoe
    x_values=np.arange(-2,-0.5,0.5)
    y_values=np.arange(0.005,0.5,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color=myColor,arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'go', markersize=1)

    # Near Wake
    x_values=np.arange(1,2,0.5)
    y_values=np.arange(0.1,4,0.2)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color=myColor,arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)

    # Far wake
    x_values=np.arange(7,8,1)
    y_values=np.arange(0.1,4,0.2)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color=myColor,arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ko',markersize=1)

    # Added to obtain saddle point for re=100
    x_values=np.arange(2,2.75,0.25)
    y_values=np.arange(1.25,2.75,0.25)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=30.0,color=myColor,arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)


    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=2,maxlength=100.0,color='black',arrowsize=0.4)

    plt.plot(1.08,3.57,color='purple',marker='P',markersize=8)
    plt.text(1.2,3.8,'D',color='purple',fontsize=16,weight='black')
    #plt.plot(0.596,0.873, color='purple',marker='P',markersize=8)
    plt.plot(1.98, 1.89, color='purple',marker='P',markersize=8)
    plt.text(2.28,2.09,'I',color='purple',fontsize=16,weight='black')
    #plt.plot(2.44, 0.36, color='purple', marker='P', markersize=8)
    #plt.plot(-0.98, 0.157, color='purple',marker='P',markersize=5)

    rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('zplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)

    plt.close()


def plot_streamlines_zplane_500(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='rainbow',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',figname='zplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-5,10))
    plt.ylim((0,6))

    colors = plt.cm.jet(np.linspace(0,1,4))
    #plt.contourf(X, Y, P, 256,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    # Upstream
    x_values=np.arange(-5,-4,0.5)
    y_values=np.arange(0.1,6,0.5)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=15,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # Horseshoe
    x_values=np.arange(-2,-0.5,0.5)
    y_values=np.arange(0.005,0.5,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=10,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'go', markersize=1)

    # Near Wake
    x_values=np.arange(1,2,0.5)
    y_values=np.arange(0.1,4,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=25,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)

    # Far wake
    x_values=np.arange(7,8,1)
    y_values=np.arange(0.1,4,0.2)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=5,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ko',markersize=1)

    # Top of cylinder
    #x_values=np.arange(0.2,0.4,0.1)
    #y_values=np.arange(4.0,4.2,0.05)
    #x,y=np.meshgrid(x_values,y_values)
    #x = x.reshape(x.size)
    #y = y.reshape(y.size)
    #seed_points = np.array([x,y])
    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=20,start_points=seed_points.T,maxlength=30.0,color='black',arrowsize=0.4)
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)

    # Added to obtain saddle point for re=100
    x_values=np.arange(2,2.75,0.25)
    y_values=np.arange(0.5,2.5,0.25)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=20,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)

    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=2,maxlength=100.0,color='black',arrowsize=0.4)
    rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)


    plt.plot(1.01,3.57,color='purple',marker='P',markersize=8)
    plt.text(1.21,3.8,'D',color='purple',fontsize=16,weight='black')
    #plt.plot(0.58,0.313, color='purple',marker='P',markersize=8)
    plt.plot(2.07, 0.88, color='purple',marker='P',markersize=8)
    plt.text(2.27,1.08,'I',color='purple',fontsize=16,weight='black')
    #plt.plot(-1.03, 0.138, color='purple',marker='P',markersize=5)
    #plt.plot(-1.46, 0.091, color='purple',marker='P',markersize=5)

    #plt.plot(3.059,0.122, color='purple',marker='P',markersize=8)
    #plt.plot(0.408, 4.097, color='purple',marker='P',markersize=8)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('zplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)

    plt.close()


def plot_streamlines_zplane_1000(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='rainbow',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',figname='zplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 12})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-5,10))
    plt.ylim((0,6))

    colors = plt.cm.jet(np.linspace(0,1,4))
    #plt.contourf(X, Y, P, 256,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    # Upstream
    x_values=np.arange(-5,-4,0.5)
    y_values=np.arange(0.1,6,0.5)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # Horseshoe
    #x_values=np.arange(-2,-0.5,0.5)
    x_values=np.full(1,-1)
    y_values=np.arange(0.005,0.5,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'go', markersize=1)

    # Near Wake
    #x_values=np.arange(1,2,0.5)
    x_values=np.full((1),1.0)
    y_values=np.arange(0.1,4,0.1)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)

    # Far wake
    x_values=np.arange(7,8,1)
    y_values=np.arange(0.3,6,0.5)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ko',markersize=1)

    # saddle point
    #1.976, 0.442
    x_values=np.full(1,1.9)
    y_values=np.arange(0.38,0.47,0.04)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=20,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4)
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)

    x_values=np.full(1,2.1)
    y_values=np.arange(0.38,0.57,0.06)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=20,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4)
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)

    # Added to obtain saddle point for re=100
    #x_values=np.arange(2,2.75,0.25)
    x_values=np.full((1),2.75)
    y_values=np.arange(0.05,2.755,0.25)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=30.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo',markersize=1)


    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=2,maxlength=100.0,color='black',arrowsize=0.4)
    rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)


    plt.plot(0.978,3.61,color='purple',marker='P',markersize=8)
    plt.text(1.2,3.8,'D',color='purple',fontsize=16,weight='black')
    #plt.annotate('D',xy=(0.978,3.61),xycoords='data',xytext=(5,5),textcoords='offset points')
    #plt.plot(0.58,0.313, color='blue', marker='X', markersize=5)
    plt.plot(1.976, 0.442, color='purple', marker='P', markersize=8)
    plt.text(2.2,0.7,'I',color='purple',fontsize=16,weight='black')
    #plt.plot(-1.103, 0.107, color='purple', marker='P', markersize=5)
    #plt.plot(-1.46, 0.091, color='black', marker='s', markersize=3)

    #plt.plot(2.522,0.079, color='purple', marker='P', markersize=8)
    #plt.plot(0.408, 4.097, color='purple', marker='P', markersize=8)

    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('zplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)

    plt.close()


def plot_streamlines_yplane_100(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='jet',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',name='4',figname='yplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-4,10))
    plt.ylim((-3,3))

    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))
    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    #upstream
    x_values=np.full(1,-2)
    y_values=np.arange(-3,3,0.4)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='black',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # near wake
    #x_values = np.arange(1.2,1.65, 0.4)
    x_values=np.full(1,1.5)
    y_values = np.arange(-0.5, 0.6, 0.1)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='black',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    x_values=np.full(1,3)
    y_values = np.arange(-0.8, 0.8, 0.1)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='black',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    x_values=np.full(1,5)
    y_values = np.arange(-0.75, 0.75, 0.1)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='black',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    plt.plot(3.06, 0.73, color='purple', marker='P', markersize=5)
    plt.plot(3.06, -0.73, color='purple', marker='P', markersize=5)

    plt.text(3.06,1.03,'A',color='purple',fontsize=16,weight='black')
    plt.text(3.06, -1.33,'B',color='purple',fontsize=16,weight='black')

    rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    #cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('yplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_streamlines_yplane_250(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='jet',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',name='4',figname='yplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-4,10))
    plt.ylim((-3,3))

    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))
    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    #upstream
    x_values=np.full(1,-2)
    y_values=np.arange(-3,3,0.4)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='black',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # near wake
    #x_values = np.arange(1.2,1.65, 0.4)
    x_values=np.full(1,1.5)
    y_values = np.arange(-0.5, 0.6, 0.25)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='black',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    plt.plot(1.41, 0.29, color='purple', marker='P', markersize=5)
    plt.plot(1.41, -0.29, color='purple', marker='P', markersize=5)

    plt.text(1.41,0.6,'A',color='purple',fontsize=16,weight='black')
    plt.text(1.41, -0.8,'B',color='purple',fontsize=16,weight='black')

    rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    #cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('yplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_streamlines_yplane_500(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='jet',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',name='4',figname='yplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-4,10))
    plt.ylim((-3,3))

    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))
    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    #upstream
    x_values=np.full(1,-2)
    y_values=np.arange(-3,3,0.4)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # near wake
    #x_values = np.arange(1.2,1.65, 0.4)
    x_values=np.full(1,1.2)
    y_values = np.arange(-0.5, 0.6, 0.2)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='dimgrey',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    # near wake
    #x_values = np.arange(1.2,1.65, 0.4)
    x_values=np.full(1,2)
    y_values = np.arange(-0.2, 0.3, 0.2)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='dimgrey',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    plt.plot(1.17, 0.33, color='purple', marker='P', markersize=5)
    plt.plot(1.17, -0.33, color='purple', marker='P', markersize=5)

    plt.text(1.17,0.7,'A',color='purple',fontsize=16,weight='black')
    plt.text(1.17, -0.9,'B',color='purple',fontsize=16,weight='black')

    rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    #cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('yplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_streamlines_yplane_1000(X,Y,U,V,P,figuresize=(7.54,3.5),colormap='jet',vmin=-0.8,vmax=0.8,step=0.2,label='Pressure',name='4',figname='yplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    plt.xlim((-4,10))
    plt.ylim((-3,3))

    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))
    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)

    #upstream
    x_values=np.full(1,-2)
    y_values=np.arange(-3,3,0.4)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])
    plt.streamplot(X, Y, U, V,linewidth=0.1,density=50,start_points=seed_points.T,maxlength=50.0,color='dimgrey',arrowsize=0.4,arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    # near wake
    #x_values = np.arange(1.2,1.65, 0.4)
    x_values=np.full(1,1.5)
    y_values = np.arange(-0.5, 0.6, 0.2)
    x, y = np.meshgrid(x_values, y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x, y])
    plt.streamplot(X, Y, U, V, linewidth=0.1, density=50, start_points=seed_points.T, maxlength=50.0, color='dimgrey',
                   arrowsize=0.4, arrowstyle='->')
    #plt.plot(seed_points[0], seed_points[1], 'ro', markersize=1)

    plt.plot(1.05, 0.36, color='purple', marker='P', markersize=5)
    plt.plot(1.05, -0.36, color='purple', marker='P', markersize=5)

    plt.text(1.17, 0.7, 'A', color='purple', fontsize=16, weight='black')
    plt.text(1.17, -0.9, 'B', color='purple', fontsize=16, weight='black')

    rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.show()
    plt.close()

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(3.5, 0.8))
    ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    #cmap = mpl.cm.jet
    cmap=plt.get_cmap(colormap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_ticks(cbarlabels, update_ticks=True)
    cb1.set_label('Pressure')
    plt.savefig('yplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)
    plt.close()


def plot_streamlines_xplane(X,Y,U,V,figuresize=(7.54,3.5),colormap='rainbow',vmin=-0.8,vmax=0.8,step=0.2,figname='yplane_stream.pdf'):

    # Set fonts
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 10})
    rc('text', usetex=True)

    fig,ax = plt.subplots(figsize=figuresize)

    cbarlabels = np.arange(vmin, vmax + step, step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i], 2)

    #plt.xlim((-4,4))
    #plt.ylim((0,5))

    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))
    #plt.contourf(X, Y, P, 100,cmap=plt.get_cmap(colormap))

    #m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    #m.set_clim(vmin, vmax)
    #m.set_array(P)
    #cbar = plt.colorbar(m)
    #cbar.set_ticks(cbarlabels)
    #cbar.set_ticklabels(cbarlabels)
    #cbar.set_label(label)



    x_values=np.arange(0.5,2,0.05)
    y_values=np.arange(-0.5,0.5,0.05)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])

    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=30,start_points=seed_points.T,maxlength=20.0,color='black',arrowsize=0.4)
    #plt.plot(seed_points[0], seed_points[1], 'ro',markersize=1)

    x_values=np.arange(0.5,6,0.2)
    y_values=np.arange(-1,1,0.05)
    x,y=np.meshgrid(x_values,y_values)
    x = x.reshape(x.size)
    y = y.reshape(y.size)
    seed_points = np.array([x,y])


    #plt.streamplot(X, Y, U, V,linewidth=0.1,density=20,start_points=seed_points.T,maxlength=20.0,color='black',arrowsize=0.4)
    #plt.plot(seed_points[0], seed_points[1], 'bo', markersize=1)

    plt.streamplot(X, Y, U, V,linewidth=0.1,density=3,maxlength=100.0,color='black',arrowsize=0.4)



    rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
    ax.add_patch(rect)
    plt.savefig(figname,bbox_inches='tight',dpi=300)
    plt.close()

    # Make a figure and axes with dimensions as desired.
    #fig = plt.figure(figsize=(3.5, 0.8))
    #ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    #cmap = mpl.cm.jet
    #cmap=plt.get_cmap(colormap)
    #norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    #cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,norm=norm,orientation='horizontal')
    #cb1.set_ticks(cbarlabels, update_ticks=True)
    #cb1.set_label('Pressure')
    #plt.savefig('yplane_stream_cbar.pdf', bbox_inches='tight', dpi=300)
    #plt.close()


def plot_average(X,Y,Z,figname,plane,xlim,ylim,bands=100,figuresize=(7.54,2.5),colormap='jet',label='Colorbar',vmin=-1,vmax=1,step=0.1):
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
    plt.contourf(X, Y, Z, bands, cmap=plt.get_cmap(colormap),vmin=v_min, vmax=v_max)
    plt.contourf(X, Y, Z, bands, cmap=plt.get_cmap(colormap), vmin=v_min, vmax=v_max)

    m = plt.cm.ScalarMappable(cmap=plt.get_cmap(colormap))
    m.set_clim(v_min, v_max)
    m.set_array(Z)
    plt.xlim(xlim)
    plt.ylim(ylim)

    cbarlabels = np.arange(v_min, v_max+step,step)
    for i in range(len(cbarlabels)):
        cbarlabels[i] = round(cbarlabels[i],2)

    cbar = plt.colorbar(m)
    cbar.set_ticks(cbarlabels)
    cbar.set_ticklabels(cbarlabels)
    cbar.set_label(label)

    #C = plt.contour(X, Y, Z, 50, linewidths=0.5 ,colors='grey',alpha=0.6)
    # C = plt.contour(X, Y, Z, 30, linewidths=0.2, colors='black')
    # plt.clabel(C, inline=1, fontsize=10)

    # Add the patch to the Axes
    if plane == 'y':
        rect = patches.Rectangle((-0.5, -0.5), 1, 1, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
        ax.add_patch(rect)
    elif plane == 'z':
        rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', facecolor='grey', alpha=0.75, zorder=10)
        ax.add_patch(rect)
    else:
        rect = patches.Rectangle((-0.5, 0), 1, 4, edgecolor='grey', alpha=0.75,linestyle='--',Fill=False)
        ax.add_patch(rect)
        pass

    plt.savefig(figname, bbox_inches='tight', dpi=300)
    plt.savefig(figname[0:-3]+'png', bbox_inches='tight', dpi=300)
    plt.close()


def write_vtk_file():
    from pyevtk.hl import gridToVTK
    from pyevtk.hl import pointsToVTK
    from pyevtk.hl import imageToVTK

    nx=1000
    ny=1000
    nz=1



    print("Writing VTK File")
    ux = np.loadtxt("Ux_zplane")
    ux = ux.reshape((nx, ny,nz))
    ux=np.dstack(ux)


    uy = np.loadtxt("Uy_zplane")
    uy = uy.reshape((nx, ny,nz))
    uy = np.dstack(uy)

    npes = get_num_proc()
    if BINARY:
        data = get_data('zplane_0', npes, prop='pmean')
    else:
        data = get_data('zplane_0.csv', npes, prop='pmean')

    X = data['X'];
    Y = data['Y'];

    xic = np.linspace(min(X), max(X), nx+1)
    yic = np.linspace(min(Y), max(Y), ny+1)
    zic = np.zeros(nz+1)

    imageToVTK("image", cellData={"ux": ux,"uy":uy},)

    x=np.zeros((nx+1,ny+1,nz+1))
    y=np.zeros((nx+1,ny+1,nz+1))
    z=np.zeros((nx+1,ny+1,nz+1))

    for k in range(nz+1):
        print(k)
        for j in range(nx+1):
            for i in range(ny+1):
                x[i,j,k]=xic[i]
                y[i,j,k]=yic[j]
                z[i,j,k]=zic[k]

    gridToVTK("grid",x,y,z,cellData={'ux':ux,'uy':uy})


def stream_fuction(x,y,u,v):
    from scipy import integrate

    # make some data
    X,Y=np.meshgrid(x,y)

    u = np.nan_to_num(u)
    v = np.nan_to_num(v)
    print(u)

    U=u
    V=v

    # integrate
    intx=integrate.cumtrapz(V,X,axis=1,initial=0)[0]
    inty=integrate.cumtrapz(U,Y,axis=0,initial=0)


    psi1=-intx+inty

    intx2=integrate.cumtrapz(V,X,axis=1,initial=0)
    inty2=integrate.cumtrapz(U,Y,axis=0,initial=0)[:,0][:,None]

    psi2=-intx2+inty2

    psi=(psi1+psi2)/2.

    u2=np.gradient(psi,axis=0)
    v2=-np.gradient(psi,axis=1)
    dx=np.gradient(X,axis=1)
    dy=np.gradient(Y,axis=0)

    u2=u2/dy
    v2=v2/dx

    plt.contourf(X, Y, psi1, 256, cmap=plt.get_cmap('RdBu'))
    plt.contour(X, Y, psi1, 20, linewidths=0.2, colors='black', alpha=0.2)
    plt.show()


if __name__ == '__main__':

    # Heading
    print('\n{:=^30}'.format('Plotting Surface'))

    # Check if parallel
    parallel = check_parallel()

    # Get Number of processors
    if (parallel):
        npes = get_num_proc()
        print("Number of Processors: {0}\n".format(npes))

    # Write VTk file
    if(False):
        write_vtk_file()

    if(False):
        if BINARY:
            data = get_data('zplane_0', npes, prop='pmean')
        else:
            data = get_data('zplane_0.csv', npes, prop='pmean')
        X = data['X'];Y = data['Y'];
        xic = np.linspace(min(X), max(X), 1000)
        yic = np.linspace(min(Y), max(Y), 1000)
        Uic = load_from_file('ux_zplane')
        Vic = load_from_file('uy_zplane')
        stream_fuction(xic,yic,Uic,Vic)

    # Instantaneous value
    if(False):
        print ('{:=^30}\n'.format('Plotting Instantaneous'))
        plane_name = "zplane"
        value = 'Ux'
        time = '100'

        inFile = plane_name + '_' + time
        outFile = plane_name + '_' + time + '_' + value + ".pdf"

        # Get the data
        data = get_data(inFile,npes,prop='Ux')
        # Plane lookup
        X_name, Y_name = plane_lookup[plane_name[0]]
        if plane_name[0] == 'x':
            X = data[Y_name];
            Y = data[X_name];
            Z = data['C']
        else:
            X = data[X_name];
            Y = data[Y_name];
            Z = data['C']

        # Create grid for interpolations
        xic = np.linspace(min(X), max(X), 1000)
        yic = np.linspace(min(Y), max(Y), 1000)

        Uic = griddata((X, Y), Z, (xic[None, :], yic[:, None]))

        # Get the data
        data = get_data(inFile, npes, prop='Uy')
        # Plane lookup
        X_name, Y_name = plane_lookup[plane_name[0]]
        if plane_name[0] == 'x':
            X = data[Y_name];
            Y = data[X_name];
            Z = data['C']
        else:
            X = data[X_name];
            Y = data[Y_name];
            Z = data['C']

        Vic = griddata((X, Y), Z, (xic[None, :], yic[:, None]))
        Pic = 0
        plot_streamlines_zplane(xic, yic, Uic, Vic, Pic, figname='zplane_stream_Ux_clean.pdf')

        # Plot the Instantaneous flow
        plot_instantaneous(data,outFile,plane=plane_name[0],bands=256,v_min=-0.8,v_max=1.5,step=0.25,c_label='ux',colormap='jet')

    # Write Values to file
    if(False):
        print('\n{:=^30}'.format('Write to File'))
        value_dict = {'wx': 'Vorticity','wy':'Vorticity, Y','Ux':'Mean Velocity','pmean':'Pressure','Uy':'Velocity, Y','Uz':'Velocity, Z'}
        #value_list=['wx','wy','ux','uy','uz','Ux','Uy','Uz','pmean','p']
        value_list=['ux','uz','uy']
        #plane_list = ['xplane1','xplane2','xplane3','xplane4','yplane1','yplane2','yplane3','yplane4','zplane']
        plane_list = ['yplane2','yplane3','zplane']

        for j in value_list:
            for i in plane_list:
                ##############################
                plane_name = i
                value=j
                clabel=j
                #######################

                # get list of time directories
                dict = get_file_list(plane_name, npes,verbose=True)

                if not os.path.isfile(value+'_'+plane_name):
                    print("\nWriting {0}".format(value+'_'+plane_name))
                    dict = create_average(dict, plane_name,prop=value)
                    write_to_file(dict['zic'],value+'_'+plane_name)

    ###################### Stream Lines ########################
    # Z plane Streamlines
    if(False):
        print('\n{:=^30}'.format('Z Plane Stream'))
        if BINARY:
            data = get_data('zplane_0', npes, prop='pmean')
        else:
            data = get_data('zplane_0.csv', npes, prop='pmean')
        X = data['X'];Y = data['Y'];
        xic = np.linspace(min(X), max(X), 1000)
        yic = np.linspace(min(Y), max(Y), 1000)
        Uic = load_from_file('ux_zplane')
        Vic = load_from_file('uy_zplane')
        Pic = 0#load_from_file('pmean_zplane')


        dir_path=os.getcwd().split('/')
        for i in dir_path:
            if i == 're_100':
                print("re_100")
                plot_streamlines_zplane_100(xic,yic,Uic,Vic,Pic,figname='zplane_stream_clean.pdf')
                break

            if i == 're_250':
                print("re_250")
                plot_streamlines_zplane_250(xic, yic, Uic, Vic, Pic, figname='zplane_stream_clean.pdf')
                break

            if i == 're_500':
                print("re_500")
                plot_streamlines_zplane_500(xic, yic, Uic, Vic, Pic, figname='zplane_stream_clean.pdf')
                break

            if i == 're_1000':
                print("re_1000")
                plot_streamlines_zplane_1000(xic, yic, Uic, Vic, Pic, figname='zplane_stream_clean.pdf')
                break

    # Y plane Streamlines
    if (False):
        num='2'
        print('\n{:=^30}'.format('Y Plane Stream'))
        if BINARY:
            data = get_data('yplane'+num+'_0', npes, prop='pmean')
        else:
            data = get_data('yplane'+num+'_0.csv', npes, prop='pmean')
        X = data['X'];
        Y = data['Z'];
        xic = np.linspace(min(X), max(X), 1000)
        yic = np.linspace(min(Y), max(Y), 1000)
        Uic = load_from_file('ux_yplane'+num)
        Vic = load_from_file('uz_yplane'+num)
        #Pic = load_from_file('p_yplane'+num)
        Pic=0

        dir_path = os.getcwd().split('/')
        for i in dir_path:
            if i == 're_100':
                print("re_100")
                plot_streamlines_yplane_100(xic, yic, Uic, Vic, Pic, name=num, figname='yplane'+num+'_stream_clean.pdf')
                break

            if i == 're_250':
                print("re_250")
                plot_streamlines_yplane_250(xic, yic, Uic, Vic, Pic, name=num, figname='yplane'+num+'_stream_clean.pdf')
                break

            if i == 're_500':
                print("re_500")
                plot_streamlines_yplane_500(xic, yic, Uic, Vic, Pic, name=num, figname='yplane'+num+'_stream_clean.pdf')
                break

            if i == 're_1000':
                print("re_1000")
                plot_streamlines_yplane_1000(xic, yic, Uic, Vic, Pic, name=num, figname='yplane'+num+'_stream_clean.pdf')
                break

    # X plane Streamlines (NOT WORKING)
    if (False):
        num='3'
        print('\n{:=^30}'.format('X Plane Stream'))
        data = get_data('xplane'+num+'_0', npes, prop='p')
        X = data['Z'];
        Y = data['Y'];
        xic = np.linspace(min(X), max(X), 1000)
        yic = np.linspace(min(Y), max(Y), 1000)
        Uic = load_from_file('uz_xplane'+num)
        Vic = load_from_file('uy_xplane'+num)
        #Pic = load_from_file('pmean_yplane'+num)
        plot_streamlines_xplane(xic, yic, Uic, Vic,figname='xplane'+num+'_stream_clean.pdf')

    ###################### X PLANE ########################
    # Average X plane Vorticity
    if (False):
        print('\n{:=^30}'.format('X Planes Vorticity'))
        value_dict = {'wx': 'Vorticity, X'}
        plane_list = ['xplane1','xplane2','xplane3','xplane4']

        v_min = -2.0
        v_max = 2.0
        step = 0.5

        for j in value_dict:
            for i in plane_list:
                print("X Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + "_average_" + value + ".pdf"

                if (BINARY):
                    csvFile = plane_name + '_' + '0'
                else:
                    csvFile = plane_name + '_' + '0.csv'


                data = get_data(csvFile, npes, prop=value)

                # Plane lookup
                X_name, Y_name = plane_lookup[plane_name[0]]
                X = data[Y_name];Y = data[X_name]

                zic = load_from_file(value+'_'+plane_name)

                # Create grid for interpolations
                xic = np.linspace(min(X), max(X), 1000)
                yic = np.linspace(min(Y), max(Y), 1000)

                zic = zic.transpose()
                plot_average(yic, xic, zic, outFile, plane_name[0], xlim, ylim,bands=100 ,label=clabel,vmin=v_min,vmax=v_max,step=step,colormap='seismic')

    # Instantaneous X plane Velocity
    if (False):
        print('\n{:=^30}'.format('Y Instantaneous Planes Velocity'))
        value_dict = {'ux': 'Velocity, X'}
        plane_list = ['xplane1']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                # print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='rainbow', v_min=-0.5, v_max=1.5, step=0.25,
                                   figuresize=(3.54, 1.75))

    # Instantaneous X plane Pressure
    if (False):
        print('\n{:=^30}'.format('Y Instantaneous Planes Pressure'))
        value_dict = {'p': 'Pressure'}
        plane_list = ['xplane1']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                # print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='rainbow', v_min=-1.2, v_max=0.8, step=0.4, figuresize=(3.54, 1.75))

    # Instantaneous X Vorticity
    if (False):
        print('\n{:=^30}'.format('X Instantaneous Vorticity'))
        value_dict = {'wy': 'Vorticity, Y', 'wx': 'Vorticity, X'}
        plane_list = ['xplane1', 'xplane2', 'xplane3', 'xplane4']
        time = '253.6'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='seismic', v_min=-10, v_max=10, step=5, figuresize=(3.54, 1.75))

    ###################### Y PLANE ########################
    # Y plane Vorticity Video
    if (False):
        print('\n{:=^30}'.format('Y Instantaneous Planes Vorticity'))
        value_dict = {'wx': 'Vorticity, X'}
        plane_list = ['yplane1']#, 'yplane2', 'yplane3', 'yplane4']
        time_list=[]
        #for i in  np.arange(21,40.5,0.5):
        for i in  np.arange(100.0,130.0,0.4):
            time_list.append(str(round(i,1)).replace(".0",''))
        count=1
        for j in value_dict:
            for i in plane_list:
                for k in time_list:
                    print("Y Plane {0}".format(i[-1]))
                    time = k
                    plane_name = i
                    value = j
                    clabel = value_dict[j]
                    xlim, ylim = get_limits(plane_name, value)
                    outFile = plane_name+ '_' + value + '_' + time + ".png"
                    csvFile=plane_name+'_'+time
                    print(csvFile)
                    data = get_data(csvFile, npes, prop=value)

                    plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,colormap='seismic',v_min=-10.0,v_max=10.0,step=4.0)
                    os.rename(outFile, "vorticity_video/img_"+str('%03d'%count)+'.png')
                    count = count + 1

    # Instantaneous Y plane Vorticity
    if (False):
        print('\n{:=^30}'.format('Y Instantaneous Planes Vorticity'))
        value_dict = ['wx']
        plane_list = ['yplane1' , 'yplane2', 'yplane3', 'yplane4']
        time='253.6'
        for j in value_dict:
            for i in plane_list:
                #print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if(BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='seismic', v_min=-5.0, v_max=5.0, step=10.0,bands=50)

    # Average Y plane Vorticity
    if (False):
        print('\n{:=^30}'.format('Y Planes Vorticity'))
        value_dict = {'wy': 'Vorticity, Y'}
        plane_list = ['yplane1', 'yplane2', 'yplane3', 'yplane4']

        v_min = -5.0
        v_max = 5.0
        step = 2.5

        for j in value_dict:
            for i in plane_list:
                print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + "_average_" + value + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + '0'
                else:
                    csvFile = plane_name + '_' + '0.csv'

                data = get_data(csvFile, npes, prop=value)

                # Plane lookup
                X_name, Y_name = plane_lookup[plane_name[0]]
                X = data[X_name];Y = data[Y_name];

                zic = load_from_file(value+'_'+plane_name)

                # Create grid for interpolations
                xic = np.linspace(min(X), max(X), 1000)
                yic = np.linspace(min(Y), max(Y), 1000)

                plot_average(xic, yic, zic, outFile, plane_name[0], (-4,15), (-4,4),label=clabel,vmin=v_min,vmax=v_max,step=step,colormap='bwr')

    # Instantaneous Y plane Velocity
    if (False):
        print('\n{:=^30}'.format('Y Instantaneous Planes Vorticity'))
        value_dict = {'ux':'Velocity, X/D','uy':'Velocity, Y/D','uz':'Velocity, Z/D'}
        plane_list = ['yplane1', 'yplane2', 'yplane3', 'yplane4']
        time = '260.8'
        if "re_100" in os.getcwd():
            time = '40'
        for j in value_dict:
            print("value: {0}".format(j))
            for i in plane_list:
                print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                #csvFile = plane_name + '_' + time + '.csv'
                csvFile = plane_name + '_' + time

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='rainbow', v_min=-0.5, v_max=1.5, step=0.25,figuresize=(3.54,1.75))

    # Instantaneous Y plane Pressure
    if (False):
        print('\n{:=^30}'.format('Y Instantaneous Planes Pressure'))
        value_dict = {'p': 'Pressure'}
        plane_list = ['yplane3']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                # print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if(BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='rainbow', v_min=-1.2, v_max=0.8, step=0.4, figuresize=(3.54, 1.75))

    # Average Y plane Velocity
    if (False):
        print('\n{:=^30}'.format('Y Planes Velocity'))
        value_dict = {'uz': 'Velocity, Z'}
        plane_list = ['yplane3']#, 'yplane2', 'yplane3', 'yplane4']

        v_min = -0.5
        v_max = 0.5
        step = 0.25

        for j in value_dict:
            for i in plane_list:
                print("Y Plane {0}".format(i[-1]))
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + "_average_" + value + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + '0'
                else:
                    csvFile = plane_name + '_' + '0.csv'

                data = get_data(csvFile, npes, prop=value)

                # Plane lookup
                X_name, Y_name = plane_lookup[plane_name[0]]
                X = data[X_name];
                Y = data[Y_name];

                zic = load_from_file(value + '_' + plane_name)

                # Create grid for interpolations
                xic = np.linspace(min(X), max(X), 1000)
                yic = np.linspace(min(Y), max(Y), 1000)

                plot_average(xic, yic, zic, outFile, plane_name[0], (-4, 15), (-4, 4), label=clabel, vmin=v_min,
                             vmax=v_max, step=step, colormap='jet',bands=60)

    ###################### Z PLANE ########################
    # Z plane Vorticity Video
    if (False):
        print('\n{:=^30}'.format('Z Instantaneous Planes Vorticity'))
        value_dict = {'wx': 'Vorticity, X'}
        plane_list = ['zplane']
        time_list=[]
        if not os.path.isdir("vorticity_video_z"):
            os.mkdir("vorticity_video_z")
        for i in  np.arange(100,130,0.5):
            time_list.append(str(i).replace(".0",''))
        count=1
        for j in value_dict:
            for i in plane_list:
                for k in time_list:
                    print("Z Plane")
                    time = k
                    plane_name = i
                    value = j
                    clabel = value_dict[j]
                    xlim, ylim = get_limits(plane_name, value)
                    outFile = plane_name+ '_' + value + '_' + time + ".png"
                    csvFile=plane_name+'_'+time + '.csv'

                    data = get_data(csvFile, npes, prop=value)

                    plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,colormap='seismic',v_min=-10.0,v_max=10.0,step=4.0,xlim=(-4,12),ylim=(0,5))
                    os.rename(outFile, "vorticity_video_z/img_"+str('%03d'%count)+'.png')
                    count = count + 1

    # Instantaneous Velocity
    if(True):
        print('\n{:=^30}'.format('Z Instantaneous Velocity'))
        value_dict = {'uz': 'Velocity, Z'}
        plane_list = ['zplane']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '112'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='jet', v_min=-1.0, v_max=1.0, step=0.5, figuresize=(3.54, 1.75))

    # Instantaneous Pressure
    if (False):
        print('\n{:=^30}'.format('Z Instantaneous Pressure'))
        value_dict = {'p': 'Pressure'}
        plane_list = ['zplane']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='rainbow', v_min=-0.8, v_max=0.8, step=0.2, figuresize=(3.54, 1.75))

    # Instantaneous Vorticity
    if (False):
        print('\n{:=^30}'.format('Z Instantaneous Vorticity'))
        value_dict = {'wy': 'Vorticity, Y','wx': 'Vorticity, X'}
        plane_list = ['zplane']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='seismic', v_min=-10, v_max=10, step=5, figuresize=(3.54, 1.75))

    # Instantaneous Velocity
    if (False):
        print('\n{:=^30}'.format('Z Instantaneous UU'))
        value_dict = {'UU': 'Velocity, X'}
        plane_list = ['zplane']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='rainbow',
                                   figuresize=(3.54, 1.75))

    # Average velocity
    if (False):
        print('\n{:=^30}'.format('Z Plane Velocity'))
        value_dict = {'Ux': 'Velocity, X'}
        plane_list = ['zplane']

        v_min = -0.6
        v_max = 1.4
        step = 0.4

        for j in value_dict:
            for i in plane_list:
                print("z Plane average")
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + "_average_" + value + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + '0'
                else:
                    csvFile = plane_name + '_' + '0.csv'

                data = get_data(csvFile, npes, prop=value)

                # Plane lookup
                X_name, Y_name = plane_lookup[plane_name[0]]
                X = data[X_name];Y = data[Y_name];

                zic = load_from_file(value+'_'+plane_name)

                # Create grid for interpolations
                xic = np.linspace(min(X), max(X), 1000)
                yic = np.linspace(min(Y), max(Y), 1000)

                plot_average(xic, yic, zic, outFile, plane_name[0], (-4,15), (0,6),label=clabel,vmin=v_min,vmax=v_max,step=step,colormap='jet')

    ###################### Turbulence ########################
    # UU
    if (False):
        print('\n{:=^30}'.format('UU Y plane'))
        value_dict = {'UU': 'UU'}
        plane_list = ['zplane']  # , 'yplane2', 'yplane3', 'yplane4']
        time = '120'
        if "re_100/" in os.getcwd():
            time = '40'
        for j in value_dict:
            for i in plane_list:
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + '_' + value + '_' + time + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + time
                else:
                    csvFile = plane_name + '_' + time + '.csv'

                data = get_data(csvFile, npes, prop=value)

                plot_instantaneous(data, outFile, plane=plane_name[0], c_label=clabel,
                                   colormap='jet', v_min=0.0, v_max=0.05, step=0.05, figuresize=(3.54, 1.75))

    ###################### Ux plots ########################
    # Average velocity
    if (False):
        print('\n{:=^30}'.format('Z Plane Velocity'))
        value_dict = {'Ux': 'Velocity, X'}
        plane_list = ['yplane1']

        for j in value_dict:
            for i in plane_list:
                print("Y Plane Plot")
                plane_name = i
                value = j
                clabel = value_dict[j]
                xlim, ylim = get_limits(plane_name, value)
                outFile = plane_name + "_average_" + value + ".pdf"
                if (BINARY):
                    csvFile = plane_name + '_' + '0'
                else:
                    csvFile = plane_name + '_' + '0.csv'

                data = get_data(csvFile, npes, prop=value)

                # Plane lookup
                X_name, Y_name = plane_lookup[plane_name[0]]
                X = data[X_name];Y = data[Y_name];

                zic = load_from_file(value + '_' + plane_name)

                # Create grid for interpolations
                xic = np.linspace(min(X), max(X), 1000)
                yic = np.linspace(min(Y), max(Y), 1000)

                Z = []

                for i in range(len(xic)):
                    if ((abs(xic[i] - 10.0)) < 1e-2):
                        column=i
                        break

                #print Z
                plt.plot(yic,zic[:,column])
                plt.show()


                xnew = np.linspace(min(Y),max(Y),90) #300 represents number of points to make between T.min and T.max

                Z=np.array(zic[1:-1,column])
                Y=np.array(yic[1:-1])


                power_smooth = spline(Y,Z,xnew,order=7)

                plt.ylim((0.6,1.1))
                plt.plot(xnew,power_smooth)
                plt.show()
