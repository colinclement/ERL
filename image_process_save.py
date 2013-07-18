#! /usr/local/enthought/bin/python

"""
Script uses image_analysis.py to load raw data, remove the background,
interpolate all images to the same scope and scale and grid. It then makes
sparse matrices of each object and saves them together in a pickled file.
"""

from image_analysis import *

super_grid = make_super_grid()

data = make_sparse_matrix(0)
M = data['sparse_data']
emittance = np.array(data['emittance'])

print 'Progress = {} out of {}'.format(0, max_number-1)

crapdata = []

for num in xrange(1,max_number):
    data = make_sparse_matrix(num)
    if data is None:
        print '****FILE NUMBER {} HAS ZERO MASS****'.format(num)
        crapdata.append(num)
        continue
    M = np.vstack((M,data['sparse_data']))
    emittance = np.hstack((emittance, data['emittance']))
    print 'Progress = {} out of {}'.format(num, max_number-1)

print 'crap data = ', crapdata

file0 = open('B1hor_processed_sparse.npy','w+')
np.savez(file0, sparse_data = M, emittances = emittance, super_grid =
         super_grid)
file0.close()
