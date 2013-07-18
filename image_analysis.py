#! /usr/local/enthought/bin/python 

"""
Script for converting .mat files into a more useful format
more to come...
Can be called inside ipython: %run matconvert (same directory)

If not called as main, super_grid = make_super_grid must be called
"""

import scipy as sp
from scipy.io import loadmat
from scipy.interpolate import griddata
from scipy.linalg import det
from scipy import sparse
import numpy as np
from numpy import sum, sqrt, einsum
import os
import matplotlib.pyplot as plt
import re

STUFFDIR = '/home/colinc/storage/ERL/ERLstuff'
walktree = os.walk(STUFFDIR,topdown=True, followlinks=False)
allfiles=[]
for folder in walktree:
    for files in folder[2]:
        allfiles.append(os.path.join(folder[0],files))
howmany=len(allfiles) #number of files accessible

#allfiles = [ os.path.join(STUFFDIR,dir,fl) for dir in os.listdir(STUFFDIR) 
#for fl in os.listdir(os.path.join(STUFFDIR,dir)) ] #(By Alex)

def sortfiles():
    """Takes allfiles and makes a dictionary with 4 different kinds of file 
    paths, sorts according to different measurements denoted by B1ver, B1hor,
    A4ver, A4hor."""
    B1ver=re.compile('B1ver')
    B1hor=re.compile('B1hor')
    A4ver=re.compile('A4ver')
    A4hor=re.compile('A4hor')
    bv=[]
    bh=[]
    av=[]
    ah=[]
    outdict={}
    for files in allfiles:
	if B1ver.search(files)!=None:
            bv.append(files)
	elif B1hor.search(files)!=None:
	    bh.append(files)
	elif A4ver.search(files)!=None:
	    av.append(files)
	elif A4hor.search(files)!=None:
	    ah.append(files)
    outdict['B1ver']=bv
    outdict['B1hor']=bh
    outdict['A4ver']=av
    outdict['A4hor']=ah
    return outdict

sortedfiles=sortfiles()

bvlen=len(sortedfiles['B1ver'])
bhlen=len(sortedfiles['B1hor'])
avlen=len(sortedfiles['A4ver'])
ahlen=len(sortedfiles['A4hor'])

def loadnum(num, kind = 'B1hor', res = 0):
    """Loads .mat file and outputs dictionary. First argument accepts one of 
    four kinds B1ver, B1hor, A4ver, A4hor as strings. Second argument
    takes ints."""
    loaded = loadmat(sortedfiles[kind][num])['data']
    #loaded is of type dtype (access with loaded.dtype)
    outdict = {}
    res = int( sqrt(len(loaded['x'][0,0])) )
    outdict['x'] = loaded['x'][0,0].reshape(res,res)
    outdict['y'] = loaded['y'][0,0].reshape(res,res)
    outdict['z'] = -loaded['z'][0,0].reshape(res,res)
    outdict['res'] = res
    #note outdict['z'] is inverted to correct data
    return outdict

def load_all_images(kind = 'B1hor'):
    """Select measurement kind and return all raw images"""
    length = len(sortedfiles[kind])
    outlist = []
    for fil in range(length):
        tp = loadmat(sortedfiles[kind][fil])['data']
	tpx = tp['x'][0,0]
	tpy = tp['y'][0,0]
	tpz = - tp['z'][0,0] #fixes inversion
	outlist.append([tpx, tpy, tpz])
    outfile = open('B1hor_raw.npy','w+')
    np.savez(outfile, data = outlist)
    outfile.close()
    return None

#load_all_images()

with open('B1hor_raw.npy','r') as infile:
    data = np.load(infile)
    B1hor_raw = data['data']
    del data

max_number = len(B1hor_raw)

def super_grid_parameters():
    """Find the size and resolution of a grid which contains all data"""
    all_x_range, all_y_range, all_pixel = [], [], []
    for image in range(B1hor_raw.shape[0]):
        all_x_range.append(B1hor_raw[image,0].ptp())
        all_y_range.append(B1hor_raw[image,1].ptp())
	all_pixel.append(int(sqrt(len(B1hor_raw[image,0]))))
    x_ranges = np.array(all_x_range)
    y_ranges = np.array(all_y_range)
    pixels = np.array(all_pixel)
    deltaXs = x_ranges/pixels
    deltaYs = y_ranges/pixels
    deltaX = np.median(deltaXs)/2. #factors chosen arbitrarily
    deltaY = np.median(deltaYs)/2.
    min_dX = min(deltaXs)
    min_dY = min(deltaYs)
    Lx = max(x_ranges)
    Ly = max(y_ranges)
    Nx = Lx/deltaX
    Ny = Ly/deltaY
    return {'Lx': Lx, 'Ly': Ly, 'deltaX': deltaXs, 'deltaY': deltaYs,
		    'x_ranges': x_ranges, 'y_ranges': y_ranges, 'pixel_no': pixels,
		    'deltaX': deltaX, 'deltaY': deltaY,
		    'Nx': Nx, 'Ny': Ny, 'min_dX': min_dX, 'min_dY': min_dY}

grid_params = super_grid_parameters()

def subtract_background(num):
    """Returns background subtracted from numpy array of B1hor image data"""
    temp = B1hor_raw[num]
    counts, bins = np.histogram(temp[2], bins = 300)
    pos = [pos for pos,value in enumerate(counts) if value == max(counts)]
    background_level = bins[pos[0]]
    temp[2] = temp[2]-background_level #subtract background
    temp[2][abs(temp[2])<100] = 0.0 #clip small numbers
    return temp

def make_super_grid():
    """For optimal resolution and range, leave custom_params = None. In order to
    create custom super_grid, specify tuple custom = (Lx, Ly, min_dX,
    min_dY)"""
    Lx = grid_params['Lx']
    Ly = grid_params['Ly']
    min_dX = grid_params['min_dX']
    min_dY = grid_params['min_dY']
    new_x = np.arange(-Lx/2., Lx/2., min_dX)
    new_y = np.arange(-Ly/2., Ly/2., min_dY)
    new_x_grid, new_y_grid = np.meshgrid(new_x, new_y)
    return np.array([new_x_grid, new_y_grid])

super_grid = make_super_grid()

def interpolate_to_center(num):
    """Select number from 0 to max_number-1. Outputs interpolated data points
    associated to super_grid specified"""
    Lx = grid_params['Lx']
    Ly = grid_params['Ly']
    min_dX = grid_params['min_dX']
    min_dY = grid_params['min_dY']
    deltaX = grid_params['deltaX']
    deltaY = grid_params['deltaY']
    indata = subtract_background(num)
    measuredpoints = np.vstack((indata[0].squeeze(),indata[1].squeeze())).T
    x_space = np.arange(indata[0].min(), indata[0].max(), deltaX)
    y_space = np.arange(indata[1].min(), indata[1].max(), deltaY)
    x_grid, y_grid = np.meshgrid(x_space, y_space)
    try:
        z_grid = griddata(measuredpoints, indata[2], (x_grid, y_grid), 
			             method = 'linear', fill_value = 0).squeeze()
    except:
	print "**First Interpolation failed**"
	raise
    mass = sum(z_grid)
    #print mass
    if mass == 0.0:
        print 'Zero mass error'
        return None
    x_center = sum(x_grid * z_grid)/mass
    y_center = sum(y_grid * z_grid)/mass
    #print x_center, y_center
    #new_x = np.arange(x_center - Lx/2., x_center + Lx/2., min_dX)
    #new_y = np.arange(y_center - Ly/2., y_center + Ly/2., min_dY)
    #new_x_grid, new_y_grid = np.meshgrid(new_x, new_y)
    new_x_grid = super_grid[0] + x_center
    new_y_grid = super_grid[1] + y_center
    try:
        z_grid = griddata(measuredpoints, indata[2], (new_x_grid, new_y_grid), 
			             method = 'linear', fill_value = 0).squeeze()
    except:
	    print "**Second Interpolation failed**"
	    raise
    #mass = sum(z_grid)
    #x_center = sum(new_x_grid * z_grid)/mass
    #y_center = sum(new_y_grid * z_grid)/mass
    #new_x_grid = new_x_grid - x_center
    #new_y_grid = new_y_grid - y_center
    return z_grid

def emittance(data):
    """Calculate the r.m.s. emittance of a uniform [X,Y,F(X,Y)] data array"""
    #tp = interpolate_to_center(num)
    xy = super_grid
    F = data 
    norm = abs(F.sum())
    return 4*sqrt(det( 1./norm * einsum('iXY,jXY,XY->ij',xy,xy,F) ))

def make_sparse_array(num):
    """Select num in range from 0 to len(B1hor_raw) = 2074 (at some time).
    Outputs sparse array and accompanying emittance in a dictionary."""
    data = interpolate_to_center(num)
    emit = emittance(data)
    if data is None:
        return None
    pos = data.nonzero()
    values = [data[pos[0][i],pos[0][i]] for i in range(len(pos[0]))]
    sparse_data = sparse.coo_matrix((values, pos), shape = data.shape)
    return {'sparse_data': sparse_data.tocsc(), 'emittance': emit}

def image_interpolated_plot(num):
    """Select number and super_grid parameters, will plot result"""
    data  = interpolate_to_center(num) 
    grid = super_grid
    if plt.fignum_exists(num): 
	    plt.close(num) 
    f1=plt.figure(num)
    ax=f1.add_subplot(111)
    cax1=ax.pcolormesh(super_grid[0], super_grid[1],data)
    f1.colorbar(cax1)
    plt.show()
    print tp[0].shape
    return f1
