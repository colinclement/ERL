#! /usr/local/enthought/bin.python

"""
Script for constructing sparse matrix of images for performing and SVD.
"""

import image_analysis as ia
import numpy as np

#Construct a list of values, rows, and columns for specifying the sparse matrix

data = ia.interpolate_to_center(0)
l, w = data.shape
arr = data.reshape(l*w)
columns = arr.nonzero()[0]
values = np.array([arr[i] for i in pos[0]])
rows = np.array([0 for i in range(len(vals))])

for num in xrange(1, 3):#ia.max_number):
    data = ia.interpolat_to_center(num)
    arr = data.reshape(l*w)
    pos = arr.nonzero()[0]
    val = np.array([arr1[i] for i in pos[0]])
    row = np.array([num for i in range(len(vals))])
    columns = np.append(columns, pos)
    rows = np.append(rows, row)
    values = np.append(values, val)

