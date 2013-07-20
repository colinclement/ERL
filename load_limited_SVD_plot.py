#! /usr/local/pub/enthought/bin/python

import numpy as np
import matplotlib.pyplot as plt
try:
    import cPickle as pickle
except ImportError:
    import pickle

with open('B1hor_limited/B1hor_limited_U_matrix.pkl','rb') as infile:
    U = pickle.load(infile)

with open('B1hor_limited/B1hor_limited_S.pkl','rb') as infile2:
    S = pickle.load(infile2)

with open('B1hor_limited/B1hor_limited_emittances.pkl','rb') as infile3:
    emittances = pickle.load(infile3)

#def image(num):
#	return V[:,num].reshape(150,150)

#def svplot(cmin = None, cmax = None, rows = 2, cols = 3):
#	fig,axes=plt.subplots(nrows = rows, ncols = cols)
#	data=np.asarray([image(n) for n in range(0, rows * cols + 1)])
#	for dat,ax in zip(data,axes.flat):
#		im=ax.imshow(dat,vmin = cmin, vmax = cmax)
#	caxes=fig.add_axes([0.91,0.1,0.02,0.8])
#	fig.colorbar(im,cax=caxes)

def projection_plot():
    """"""
    figs, axs = plt.subplots(nrows = 4, ncols = 4, figsize=(18,18))
    colors = np.log(emittances+0.005)
    cm = plt.cm.get_cmap('RdBu') #PiYg
    projs = U #* S[:,None]
    for row in range(4):
        for col in range(4):
            if col>=row:
                sc = axs[row,col].scatter(projs[:,row], projs[:,col+1], alpha = 0.55,
		c = colors, marker = 'x', cmap = cm, vmin = min(colors), vmax = max(colors));
                axs[row,col].set_title('{} versus {}'.format(row+1,col+2))
            else:
                axs[row,col].axis('off');
    caxes = figs.add_axes([0.4, 0.3, 0.02, 0.19])
    plt.colorbar(sc, cax = caxes)
    return None

