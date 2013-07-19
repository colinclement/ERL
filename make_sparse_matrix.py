#! /usr/local/enthought/bin.python

"""
Script for constructing sparse matrix of images for performing and SVD.
"""

import image_analysis as ia
import numpy as np
import scipy.sparse as sparse
try:
    import cPickle as pickle
except ImportError:
    import pickle

#Construct a list of values, rows, and columns for specifying the sparse matrix

data = ia.interpolate_to_center(0)
emittances = np.array(ia.emittance(data))
l, w = data.shape
arr = data.reshape(l*w)
columns, = arr.nonzero()
values = np.array([arr[i] for i in columns])
rows = np.array([0 for i in range(len(values))])
print 'Progress 0 out of {}'.format(ia.max_number-1)

for num in xrange(1, ia.max_number):
    data = ia.interpolate_to_center(num)
    emit = np.array(ia.emittance(data))
    emittances = np.hstack((emittances, emit))
    arr = data.reshape(l*w)
    pos, = arr.nonzero()
    vals = np.array([arr[i] for i in pos])
    new_row = np.array([num for i in range(len(vals))])
    columns = np.append(columns, pos)
    rows = np.append(rows, new_row)
    values = np.append(values, vals)
    print 'Progress = {} out of {}'.format(num, ia.max_number-1)

sparse_M = sparse.coo_matrix((values, (rows, columns)), shape =
                             (4, l*w)).tocsc() #ia.max_number
outdict = {'sparse_data': sparse_M, 'emittances': emittances}

with open('B1hor_sparse.pkl','w+') as outfile:
    pickle.dump(outdict, outfile, protocol = -1)

