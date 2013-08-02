#! /usr/local/pub/enthought/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col

with open('B1hor_limited_results/B1hor_sparse_limited_SVD_p.npy','rb') as infile3:
    data = np.load(infile3)
    del data.f
    Uk = data['Uk']
    Sk = data['Sk']
    Vk = data['Vtk'].T
    emittances = data['emittances']
    l, w = data['shape']
    del data

nanvalues = [i for i in range(len(emittances)) if np.isnan(emittances[i])]     
#set nans to finite value (near minimum)
#Only 3 nan values exist so this isn't important
emittances[nanvalues] = 0.0004

def image(num):
	return Vk[:,-num].reshape(l,w)

def svplot(cmin = None, cmax = None, window=300):
    fig,axes=plt.subplots(nrows = 2, ncols = 4)
    cl=int(l/2.)
    cw=int(w/2.)
    data=np.asarray([ Vk[:,n].reshape(l,w).T[cw-window:cw+window,
                                             cl-window:cl+window] for n in range(10)[::-1] ])
    for dat,ax in zip(data,axes.flat):
        im=ax.imshow(dat,vmin = cmin, vmax = cmax)
    caxes=fig.add_axes([0.91,0.1,0.02,0.8])
    fig.colorbar(im,cax=caxes)

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

def projection_plot(alpha = 0.4, plotorder = 'lowtohigh'):
    """Creates an upper triangular array of plots of the projection of the data
    onto planes spanned by the singular vectors"""
    figs, axs = plt.subplots(nrows = 4, ncols = 4, figsize=(18,18))
    plotorderdict = {'lowtohigh': emittances.argsort()[::-1], 'hightolow':
                     emittances.argsort()}
    plotorder = plotorderdict[plotorder]
    projs = Uk.T[::-1].T #* Sk
    for row in range(4):
        for col in range(4):
            if col>=row:
                px = projs[:,row][plotorder]
                py = projs[:,col+1][plotorder]
                sc = axs[row,col].scatter(px, py,
                                          alpha = alpha,
		                                  c = emittances[plotorder], marker =
                                          'o', linewidth = 0, 
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

