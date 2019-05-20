"""
Plot NH LORRI results.
"""
import numpy as np
import matplotlib.pyplot as plt
#import png
#import itertools
from scipy.io import loadmat
#from scipy.ndimage import filters,zoom
#import skimage
#from skimage import io,color
#from PIL import Image
#import pyfits
from matplotlib.backends.backend_pdf import PdfPages

fig=plt.figure(figsize=(6,5))

ax = fig.add_subplot(1,1,1)

lorri = loadmat('../../scratch/nh_make_results.mat')

lorri['distance'][0] = lorri['distance'][0] - 0.07
lorri['distance'][1] = lorri['distance'][1] + 0.07

ax.plot(lorri['distance'],lorri['rawmean'],linestyle='',marker='o')
ax.errorbar(lorri['distance'],lorri['rawmean'],yerr=lorri['rawerr'],linestyle='',marker='o')

ax.set_ylim([0,120])
ax.set_xlim([0,20])
ax.set_xlabel('Heliocentric Distance (AU)')
ax.set_ylabel(r'$\lambda I_{\lambda}^{\rm diffuse}$ in Masked Image (nW m$^{-2}$ sr$^{-1}$)')


#plt.tight_layout()
#plt.show
pdf = PdfPages('nh_plot_raw.pdf')
pdf.savefig()
pdf.close()

fig=plt.figure(figsize=(6,5))

ax = fig.add_subplot(1,1,1)

ax.plot(lorri['distance'],lorri['cobmean'],linestyle='',marker='o',color='blue')
ax.errorbar(lorri['distance'],lorri['cobmean'],yerr=lorri['coberr'],linestyle='',marker='o',color='blue')

ax.plot([0,20],[0,0],linestyle=':',color='black')

ax.plot([0,20],[np.squeeze(lorri['supermean']),np.squeeze(lorri['supermean'])],color='blue',linestyle='--')
ax.plot([0,20],[np.squeeze(lorri['supermean']+lorri['supererr']),np.squeeze(lorri['supermean']+lorri['supererr'])],color='blue',linestyle=':')
ax.plot([0,20],[np.squeeze(lorri['supermean']-lorri['supererr']),np.squeeze(lorri['supermean']-lorri['supererr'])],color='blue',linestyle=':')

ax.set_ylim([-20,55])
ax.set_xlim([0,20])
ax.set_xlabel('Heliocentric Distance (AU)')
ax.set_ylabel(r'LORRI $\lambda I_{\lambda}^{\rm resid}$ (nW m$^{-2}$ sr$^{-1}$)')


#plt.tight_layout()
#plt.show
pdf = PdfPages('nh_plot_results.pdf')
pdf.savefig()
pdf.close()
