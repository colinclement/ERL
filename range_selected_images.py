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
with open('range_limited/B1hor_lr_super_grid.npy', 'w+') as sf:
    np.savez(sf, super_grid = ia.super_grid)
          
data = ia.interpolate_to_center(share_limit[0])
emit = np.array(ia.emittance(data))
l,w = data.shape
data = data.reshape(l*w)
print 'Progress = {} out of {}'.format(0, number_limit-1)
with open('range_limited/B1hor_lr_{}.npy'.format(share_limit[0]),'w+') as df:
    np.savez(df, M = data, emittance = emit, shape = (l,w)) 

crapdata = []
for num in share_limit[1:]:
    data = ia.interpolate_to_center(num)
    emit = ia.emittance(data)
    if data is None:
        print '****FILE NUMBER {} HAS ZERO MASS****'.format(num)
        crapdata.append(num)
        continue
    #M = np.vstack((M,data.reshape(l*w)))
    #emittance = np.hstack((emittance, np.array(emit)))
    with open('range_limited/B1hor_lr_{}.npy'.format(num),'w+') as df:
        np.savez(df, M = data, emittance = emit, shape = (l,w)) 
    print 'Progress = {} out of {}'.format(num, share_limit.max())

if crapdata != []:
    print 'crap data = ', crapdata
"""
U, S, Vt = svd(M, full_matrices = False)

file0 = open('B1hor_range_limited_processed.npy','w+')
np.tofile(file0,{'data': M, 'emittances': emittance, 'super_grid': super_grid})
file0.close()

file1 = open('B1hor_range_limit_SVD.npy','w')
np.tofile(file1,{'U': U, 'S': S,'Vt': Vt, 'emittances': emittance})
file1.close()"""
