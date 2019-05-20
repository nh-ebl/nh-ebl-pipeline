
# coding: utf-8

# In[2]:

get_ipython().magic(u'matplotlib inline')


# In[3]:

# PART 1
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

import scipy.io


# In[4]:

#load .mat file
mat = scipy.io.loadmat('nh_lorri_psf_updated.mat',matlab_compatible=True)
#print('Variable names:')
#print(mat.keys()) #print out keywords (names) of the variables

#the useful data is in key='psf'
data = mat['psf']

#in variable 'psf', we have two fields
#'psf' and 'resizepsf'
#each reside at the [0,0] location
psf = data['psf'][0,0] #native resolution
modelpsf = data['modelpsf'][0,0] #psf used for modeling
onexpsf = data['onexpsf'][0,0] #native resolution again (1x)
fourxpsf = data['fourxpsf'][0,0] #4x, native resolution for the fine sampling, used for image
tenxpsf = data['tenxpsf'][0,0] #10x
fourtyxpsf = data['fourtyxpsf'][0,0] #40x, used for radial curve
pixelwidths = data['pixelwidths'][0,0] #pixel widths array
centerpix = data['centerpix'][0,0] #position of central pixels (locations of maxima) in each psf

#extract pixel width information from pixelwidths array
pixelwidth_native = pixelwidths[0,0]
pixelwidth_model = pixelwidths[0,1]
pixelwidth_one = pixelwidths[0,2]
pixelwidth_four = pixelwidths[0,3]
pixelwidth_ten = pixelwidths[0,4]
pixelwidth_fourty = pixelwidths[0,5]

print pixelwidths
print pixelwidth_four


# In[5]:

#calculate the annular average for the fine image
x = np.around(np.arange(-13.25,13.25+0.053,0.106),decimals=3)#use the center of each pixel as the coordinate
y = np.around(np.arange(-13.25,13.25+0.053,0.106),decimals=3)
xx,yy = np.meshgrid(x[:-1],y[:-1])
r = np.around((xx**2.0 + yy**2.0)**(1.0/2),decimals=3)

r = r.reshape((62500,1))
resizepsf_1D = fourtyxpsf.reshape((62500,1))

hist, bin_edge = np.histogram(r,bins=125)
mid_r = np.zeros(len(bin_edge)-1)
average_B = np.zeros(len(bin_edge)-1)
std_B = np.zeros(len(bin_edge)-1)

for i in range(0,(len(bin_edge)-1)):
    r_low = bin_edge[i]
    r_high = bin_edge[i+1]
    outer_r_index = (r[:] < r_high).astype(int)
    inner_r_index = (r[:] >= r_low).astype(int)
    nan_index = (~np.isnan(resizepsf_1D)).astype(int)
    annulus_index = outer_r_index & inner_r_index & nan_index
    
    B_in_annulus = resizepsf_1D[annulus_index>0]
    mid_r[i] = (r_low+r_high)/2.0
    #mid_r[i] = r_low
    average_B[i] = np.average(B_in_annulus)
    std_B[i] = np.std(B_in_annulus)
    
average_B = average_B/np.max(average_B)


# In[6]:

#calculate the annular average for the coarse image
x_coarse = np.around(np.arange(-12.72,12.72+2*1.06,1.06),decimals=3) 
y_coarse = np.around(np.arange(-12.72,12.72+2*1.06,1.06),decimals=3)
xx_coarse,yy_coarse = np.meshgrid(x_coarse[:-1],y_coarse[:-1])
r_coarse = np.around((xx_coarse**2.0 + yy_coarse**2.0)**(1.0/2),decimals=3)
print x_coarse
r_coarse = r_coarse.reshape((625,1))
psf_1D = fourxpsf.reshape((625,1))

#we want the bin edge to coincide with the 1/2 pixel width
bin_edge_coarse = np.asarray([0,1,3,5,7,9,11,13,15,17,19,21,23,25])
bin_edge_coarse = bin_edge_coarse*1.06/2.0

print bin_edge_coarse
mid_r_coarse = np.zeros(len(bin_edge_coarse)-1)
average_B_coarse = np.zeros(len(bin_edge_coarse)-1)
std_B_coarse = np.zeros(len(bin_edge_coarse)-1)

for i in range(0,(len(bin_edge_coarse)-1)):
    r_low_coarse = bin_edge_coarse[i]
    r_high_coarse = bin_edge_coarse[i+1]
    outer_r_index_coarse = (r_coarse[:] < r_high_coarse).astype(int)
    inner_r_index_coarse = (r_coarse[:] >= r_low_coarse).astype(int)
    nan_index_coarse = (~np.isnan(psf_1D)).astype(int)
    annulus_index_coarse = outer_r_index_coarse & inner_r_index_coarse & nan_index_coarse
    
    B_in_annulus_coarse = psf_1D[annulus_index_coarse>0]
    mid_r_coarse[i] = (r_low_coarse+r_high_coarse)/2.0
    average_B_coarse[i] = np.average(B_in_annulus_coarse)
    std_B_coarse[i] = np.std(B_in_annulus_coarse)
    
average_B_coarse = average_B_coarse/np.max(average_B_coarse)
print bin_edge_coarse
print(average_B_coarse)


# In[86]:

#make a mesh grid to plot 2D
dx_4 = 1.06
dy_4 = dx_4
x_4,y_4 = np.mgrid[slice(-13.25,13.25+dx_4,dx_4),slice(-13.25,13.25+dy_4,dy_4)]
z_4 = np.log10(fourxpsf)

#define color scale
levels = MaxNLocator(nbins=100).tick_values(-4,0)
cmap = plt.get_cmap('gnuplot2')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#1st figure: pixel = 0.43"
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,3.))
#fig.subplots_adjust(left=-0.5)

ax1.set(adjustable='box', aspect='equal')
#ax1.set_aspect('equal')
im = ax1.pcolormesh(x_4,y_4,z_4,cmap=cmap,norm=norm)
#ax1.set_title('PSF')
ax1.set_xlabel('Arcsec')
ax1.set_ylabel('Arcsec',labelpad=-0.75)
ax1.set_xlim(-13.25,13.25)
ax1.set_ylim(-13.25,13.25)
#add a color bar
cax = fig.add_axes([0.41,0.25,0.02,0.65]) #[left margin, bottom margin, width, height]
cbar = fig.colorbar(im,ax=ax1,cax=cax,orientation='vertical',ticks=[-4,-3,-2,-1,0])
cbar.ax.set_xlabel('      log$_{10}$B',labelpad=5)


#2nd figure: annular average
#coarse image:
width = np.zeros(len(average_B_coarse))
width = width+1.06
width[0] = width[0]/2.0
ax2.bar(bin_edge_coarse[:-1],average_B_coarse,width,label='PSF, stacked image',color='#4eb3d3',alpha=0.25,zorder=0)

ax2.axvline(0.75,color='#253494',linestyle=':',linewidth=1,label='HWHM, Cheng et al 2008')
ax2.axvline(0.67/2,color='k',linestyle='--',linewidth=1,label='Diffraction Limit')
ax2.axhline(0.5,color='#d95f0e',linestyle='-')
ax2.text(6.75,0.65,'HWHM, this paper')
ax2.axvline(pixelwidth_four/2,color='#c51b8a',linestyle='--',linewidth=1.5,label='Half Pixel Width')

#fine image
ax2.fill_between(mid_r,average_B+std_B, average_B,color='#78c679',alpha=0.5)
ax2.plot(mid_r,average_B,color='#238443',linewidth=2,label='PSF, fine image',zorder=1,alpha=0.75)
ax2.fill_between(mid_r,average_B, average_B-std_B,color='#78c679',alpha=0.5)


ax2.set_xlabel('Radius (arcsec)')
ax2.set_ylabel('B($\lambda$)',labelpad=4)
ax2.set_yscale('log')
ax2.set_xlim(0.1,15)
ax2.set_ylim(0.0001,2)
ax2.legend(loc=(0.31,0.3),frameon=False,prop={'size':10},numpoints=1)

#save figure
plt.tight_layout(w_pad=3.0)
plt.savefig('psf_annular_average.pdf',bbox_inches='tight')
plt.show()

