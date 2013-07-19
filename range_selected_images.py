#! /usr/local/enthought/bin/python

import image_analysis as ia
from scipy.linalg import svd
import numpy as np

x_ranges = ia.grid_params['x_ranges']
y_ranges = ia.grid_params['y_ranges']

x_limit = 0.3
y_limit = 2.

x_below_limit = [i for i in range(len(x_ranges)) if x_ranges[i] < x_limit]
y_below_limit = [i for i in range(len(y_ranges)) if y_ranges[i] < y_limit]

share_limit = np.intersect1d(x_below_limit, y_below_limit)
number_limit = len(share_limit)

#I estimate the memory usage for the full matric with these parameters 
#would be 11GB
#change grid parameters
ia.grid_params['Lx'] = x_limit
ia.grid_params['Ly'] = y_limit
#Remake super grid
ia.super_grid = ia.make_super_grid()

data = ia.interpolate_to_center(0)
l,w = data.shape
M = data.reshape(l*w)
emittance = np.array(ia.emittance(data))
print 'Progress = {} out of {}'.format(0, number_limit-1)

crapdata = []

for num in share_limit:
    data = ia.interpolate_to_center(num)
    emit = ia.emittance(data)
    if data is None:
        print '****FILE NUMBER {} HAS ZERO MASS****'.format(num)
        crapdata.append(num)
        continue
    M = np.vstack((M,data.reshape(l*w)))
    emittance = np.hstack((emittance, np.array(emit)))
    print 'Progress = {} out of {}'.format(num, number_limit-1)

if crapdata != []:
    print 'crap data = ', crapdata

U, S, Vt = svd(M, full_matrices = False)

file0 = open('B1hor_range_limited_processed.npy','w+')
np.tofile(file0,{'data': M, 'emittances': emittance, 'super_grid': super_grid})
file0.close()

file1 = open('B1hor_range_limit_SVD.npy','w')
np.tofile(file1,{'U': U, 'S': S,'Vt': Vt, 'emittances': emittance})
file1.close()
