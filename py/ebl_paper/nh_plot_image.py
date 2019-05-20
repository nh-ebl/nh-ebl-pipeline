"""
Plot the calibrated images of each field.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.colorbar as clbr
import matplotlib.patheffects as path_effects
#import png
#import itertools
from scipy.io import loadmat
#from scipy.ndimage import filters,zoom
#import skimage
#from skimage import io,color
#from PIL import Image
#import pyfits
from matplotlib.backends.backend_pdf import PdfPages

fig=plt.figure(figsize=(6.5,5))

extent = [-4.3 * 256 / 60 / 2,4.3 * 256 / 60 / 2,\
          -4.3 * 256 / 60 / 2,4.3 * 256 / 60 / 2]

myimage = loadmat('../../scratch/field5_image58.mat')

ax1 = fig.add_subplot(2,2,1)

ax1.imshow(myimage['image'],extent=extent,clim=[-100,400])

ax1a=ax1.annotate('Field 1',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
             horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax1a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])


ax1.set_ylabel('arcmin')
ax1.set_xlabel('arcmin')

myimage = loadmat('../../scratch/field3_image46.mat')

ax2 = fig.add_subplot(2,2,2)

ax2.imshow(myimage['image'],extent=extent,clim=[-100,400])

ax2a=ax2.annotate('Field 2',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
            horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax2a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])

ax2.set_ylabel('arcmin')
ax2.set_xlabel('arcmin')

myimage = loadmat('../../scratch/field6_image61.mat')

ax3 = fig.add_subplot(2,2,3)

ax3.imshow(myimage['image'],extent=extent,clim=[-100,400])

ax3a=ax3.annotate('Field 3',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
            horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax3a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])

ax3.set_ylabel('arcmin')
ax3.set_xlabel('arcmin')

myimage = loadmat('../../scratch/field7_image64.mat')

ax4 = fig.add_subplot(2,2,4)

ax4.imshow(myimage['image'],extent=extent,clim=[-100,400])

ax4a=ax4.annotate('Field 4',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
            horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax4a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])

ax4.set_ylabel('arcmin')
ax4.set_xlabel('arcmin')

norm = clr.Normalize(vmin=-100, vmax=400)
cx1 = fig.add_axes([0.81, 0.125, 0.05, 0.84])
cb1 = clbr.ColorbarBase(cx1, #cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_label('$\lambda I_{\lambda}$ (nW m$^{-2}$ sr$^{-1}$)')

plt.tight_layout(rect=[0,0,0.85,1])
plt.show()

#pdf = PdfPages('nh_plot_image.pdf')
#pdf.savefig()
fig.savefig('nh_plot_image.pdf')
plt.close

##########################################################################
##########################################################################

fig=plt.figure(figsize=(6.5,5))

extent = [-4.3 * 256 / 60 / 2,4.3 * 256 / 60 / 2,\
          -4.3 * 256 / 60 / 2,4.3 * 256 / 60 / 2]

myimage = loadmat('../../scratch/field5_masked58.mat')
whpl = myimage['nanimage'] > 200
myimage['nanimage'][whpl] = np.nan

ax1 = fig.add_subplot(2,2,1)

ax1.imshow(myimage['nanimage'],extent=extent,clim=[-100,400])

ax1a=ax1.annotate('Field 1',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
             horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax1a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])

ax1.set_ylabel('arcmin')
ax1.set_xlabel('arcmin')

myimage = loadmat('../../scratch/field3_masked46.mat')
whpl = myimage['nanimage'] > 200
myimage['nanimage'][whpl] = np.nan

ax2 = fig.add_subplot(2,2,2)

ax2.imshow(myimage['nanimage'],extent=extent,clim=[-100,400])

ax2a=ax2.annotate('Field 2',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
            horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax2a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])

ax2.set_ylabel('arcmin')
ax2.set_xlabel('arcmin')

myimage = loadmat('../../scratch/field6_masked61.mat')
whpl = myimage['nanimage'] > 200
myimage['nanimage'][whpl] = np.nan

ax3 = fig.add_subplot(2,2,3)

ax3.imshow(myimage['nanimage'],extent=extent,clim=[-100,400])

ax3a=ax3.annotate('Field 3',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
            horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax3a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])

ax3.set_ylabel('arcmin')
ax3.set_xlabel('arcmin')

myimage = loadmat('../../scratch/field7_masked64.mat')
whpl = myimage['nanimage'] > 200
myimage['nanimage'][whpl] = np.nan

ax4 = fig.add_subplot(2,2,4)

ax4.imshow(myimage['nanimage'],extent=extent,clim=[-100,400])

ax4a=ax4.annotate('Field 4',xy=(-1,8),xycoords='data',
            xytext=(-1,8), textcoords='data',
            horizontalalignment='right', verticalalignment='top',\
                  color='yellow',fontsize=14)
ax4a.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'),
                       path_effects.Normal()])


ax4.set_ylabel('arcmin')
ax4.set_xlabel('arcmin')

norm = clr.Normalize(vmin=-100, vmax=400)
cx1 = fig.add_axes([0.81, 0.125, 0.05, 0.84])
cb1 = clbr.ColorbarBase(cx1, #cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_label('$\lambda I_{\lambda}$ (nW m$^{-2}$ sr$^{-1}$)')

plt.tight_layout(rect=[0,0,0.85,1])
plt.show()

#pdf = PdfPages('nh_plot_image.pdf')
#pdf.savefig()
fig.savefig('nh_plot_masked.pdf')
plt.close
