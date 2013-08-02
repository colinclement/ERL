#! /usr/local/enthought/bin.python

"""
Script for constructing sparse matrix of images for performing and SVD.
"""

import image_analysis as ia
import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import svds
from multiprocessing import Pool
from datetime import datetime

start_time = datetime.now()

#Construct a list of values, rows, and columns for specifying the sparse matrix
#ia.max_number = 15 #for testing purposes


def dense_values(num):
    data = ia.interpolate_to_center(num)
    emit = np.array(ia.emittance(data))
    if emit == None:
        return None
    l,w = data.shape
    arr = data.reshape(l*w)
    cols, = arr.nonzero()
    vals = arr[cols] #np.array([arr[i] for i in cols])
    row = num*np.ones((len(vals))) #np.array([num for i in range(len(vals))])
    print 'Completed = {} out of {}'.format(num, ia.max_number-1)
    return cols, vals, row, emit

if __name__ == '__main__':
    data = ia.interpolate_to_center(0)
    l,w = data.shape
    del data

    columns = np.array([])
    values = np.array([])
    rows = np.array([])
    emittances = np.array([])

    pool = Pool(processes = 16)
    results = pool.map(dense_values, range(ia.max_number))

    print 'Sparse values calculated'

    for i in range(ia.max_number):
        if results[i] == None:
            continue
        columns = np.hstack((columns, results[i][0]))
        values = np.hstack((values, results[i][1]))
        rows = np.hstack((rows, results[i][2]))
        emittances = np.hstack((emittances, results[i][3]))

    print 'Sparse values compiled'


    sparse_M = sparse.coo_matrix((values, (rows, columns)), shape =
                             (ia.max_number, l*w)).tocsc() 

    print 'Sparse matrix made'

    Uk, Sk, Vtk = svds(sparse_M, k = 6)

    

    print 'SVD complete'

    #outdict = {'Uk': Uk, 'Sk': Sk, 'Vtk': Vtk, 'emittances': emittances}

    with open('B1hor_sparse_SVD_pnormed.npy','wb') as outfile:
        #pickle.dump(outdict, outfile, protocol = -1)
        np.savez(outfile, Uk = Uk, Sk = Sk, Vtk = Vtk, emittances = emittances,
                 shape = (l,w))
    print 'Time to complete {}'.format(datetime.now()-start_time)
