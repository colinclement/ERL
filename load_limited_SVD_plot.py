#! /usr/local/pub/enthought/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
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

nanvalues = [i for i in range(len(emittances)) if np.isnan(emittances[i])]     
#set nans to finite value (near minimum)
#Only 3 nan values exist so this isn't important
emittances[nanvalues] = 0.0004

#def image(num):
#	return V[:,num].reshape(150,150)

#def svplot(cmin = None, cmax = None, rows = 2, cols = 3):
#	fig,axes=plt.subplots(nrows = rows, ncols = cols)
#	data=np.asarray([image(n) for n in range(0, rows * cols + 1)])
#	for dat,ax in zip(data,axes.flat):
#		im=ax.imshow(dat,vmin = cmin, vmax = cmax)
#	caxes=fig.add_axes([0.91,0.1,0.02,0.8])
#	fig.colorbar(im,cax=caxes)

my_cdict = {'red': ((0.0, 0.0, 0.0),
                    (0.06, 0.0, 0.0),
                    (0.1, 0.7, 0.7),
                    (0.15, 0.7, 0.7),
                    (1.0, 1.0, 1.0)),
           'green': ((0.0, 1.0, 1.0),
                     (0.06, 0.5, 0.5), 
                     (0.1, 0.0, 0.0),
                     (0.15, 0.0, 0.0),
                     (1.0, 0.0, 0.0)),
           'blue': ((0.0, 0.0, 0.0),
                    (0.06, 0.0, 0.0),
                    (0.1, 0.7, 0.7),
                    (0.15, 0.5, 0.5),
                    (1.0, 0.0, 0.0))}
my_cmap = col.LinearSegmentedColormap('my_cmap', my_cdict, N=256, gamma = 1.)

def projection_plot(alpha = 0.6):
    """Creates an upper triangular array of plots of the projection of the data
    onto planes spanned by the singular vectors"""
    figs, axs = plt.subplots(nrows = 4, ncols = 4, figsize=(18,18))
    plotorder = emittances.argsort()[::-1]
    projs = U #* S[:,None]
    for row in range(4):
        for col in range(4):
            if col>=row:
                px = projs[:,row][plotorder]
                py = projs[:,col+1][plotorder]
                sc = axs[row,col].scatter(px, py,
                                          alpha = alpha,
		                                  c = emittances[plotorder], marker = 'o', 
                                          cmap = my_cmap, vmin = emittances.min(), 
                                          vmax = emittances.max());
                axs[row,col].set_title('{} versus {}'.format(row+1,col+2))
            else:
                axs[row,col].axis('off');
    caxes = figs.add_axes([0.4, 0.15, 0.02, 0.19])
    plt.colorbar(sc, cax = caxes)
    emit_spec_axs = figs.add_axes([0.15, 0.15, 0.18, 0.18])
    emit_spec = np.sort(emittances)
    emit_spec_axs.scatter(range(len(emit_spec)),emit_spec,
                                c = emit_spec, cmap = my_cmap, linewidth=0)
    emit_spec_axs.set_title('Spectrum of emittances')
    return None

