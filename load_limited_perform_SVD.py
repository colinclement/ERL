#! /usr/local/enthought/bin/python 

"""
This script imports the range-limited images which have been processed. It
composes a matrix from them, performs an SVD, then saves the SVD to disk.
"""

import image_analysis as ia
import numpy as np
from scipy.linalg import svd
try:
    import cPickle as pickle
except ImportError:
    import pickle
import os

walktree = os.walk('range_limited/', topdown = True)
file_names = walktree.next()
file_list = [os.path.join('range_limited/',file) for file in file_names[2]]
file_list = [file for file in file_list if not file
             =='range_limited/B1hor_lr_super_grid.npy']
with open(file_list[0],'r') as infile:
    data = np.load(infile)
    l, w  = data['shape']
    M = data['M'].reshape(l*w)
    emittance = data['emittance']
    del data
index = 1
for file in file_list[1:]:
    with open(file,'r') as infile:
        data = np.load(infile)
        new_row = data['M'].reshape(l*w)
        emit = data['emittance']
        del data
    M = np.vstack((M, new_row))
    emittance = np.hstack((emittance,emit))
    index = index + 1
    print 'Loaded {} successfull'.format(file)
U, S, Vt =  svd(M, full_matrices = False, overwrite_a = True)

with open('B1hor_limited_U_matrix.pkl','w+') as outfile:
    pickle.dump(U, outfile)
    del U
with open('B1hor_limited_Vt_matrix.pkl', 'w+') as outfile:
    pickle.dump(Vt, outfile)
    del Vt
with open('B1hor_limited_S.pkl', 'w+') as outfile:
    pickle.dump(S, outfile)
    del S
with open('B1hor_limited_emittances.pkl','w+') as outfile:
    pickle.dump(emittances, outfile)

