#! /usr/local/pub/enthought/bin/python

import numpy as np
import matplotlib.pyplot as plt

with open('B1horSVD_refined1.npy','r') as infile:
    data = np.load(infile)
    U=data['U']
    S=data['S']
    V=data['Vt'].T
    del data

with open('B1horPCA.npy','r') as infile2:
    data = np.load(infile2)
    PCA_U = data['PCA_U']
    PCA_S = data['PCA_S']
    PCA_Vt = data['PCA_Vt']
    del data

with open('emittances.npy','r') as infile3:
    data = np.load(infile3)
    emittances = data['B1hor_emittances']
    del data

def image(num):
	return V[:,num].reshape(150,150)

def svplot(cmin = None, cmax = None, rows = 2, cols = 3):
	fig,axes=plt.subplots(nrows = rows, ncols = cols)
	data=np.asarray([image(n) for n in range(0, rows * cols + 1)])
	for dat,ax in zip(data,axes.flat):
		im=ax.imshow(dat,vmin = cmin, vmax = cmax)
	caxes=fig.add_axes([0.91,0.1,0.02,0.8])
	fig.colorbar(im,cax=caxes)

def data_slices_plot(PCA = True):
    """"""
    figs, axs = plt.subplots(nrows = 4, ncols = 4, figsize=(18,18))
    colors = np.log(emittances+0.005)
    cm = plt.cm.get_cmap('RdBu') #PiYg
    if PCA:
        projs = PCA_U #* PCA_S[:,None]
    else:
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
    
