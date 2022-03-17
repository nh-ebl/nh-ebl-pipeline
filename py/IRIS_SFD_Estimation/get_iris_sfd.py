################################################################################
# NAME : model_test.py
# DATE STARTED : July 5, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : The purpose of this script is to test that the map we recreate using
# the best fit params from Planck agree with their model.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
import os
from astropy.io import fits
sys.path.append('../')
from utilities import *
from make_map import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from astropy import units as u
from astropy.constants import c

# Define custom font sizes for figures
SMALL_SIZE = 21
MEDIUM_SIZE = 25
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fieldnum = 1
lam = 100*u.micron
nu = (c.to(u.micron/u.s)/lam).to(u.Hz)

#this is a reference image
# hdul = fits.open('regist_pluto_20151102_030873_lor_0308731028_0x633_sci.fit')
#<read in the fits image name from a text file made by matlab>
currentDir = os.path.dirname(os.path.realpath(__file__)) #get current working directory (calling from cmd line python has directory issues)
os.chdir(currentDir) #change dir to the needed directory
fileID = open("imagefile.txt","r") #open the file
imagefile = fileID.readlines()[0] #get the imagefile path
fileID.close() #close the file
hdul = fits.open(imagefile) #use the imagefile path

ref_head = hdul[0].header
timestamp = ref_head['SPCUTCJD'][3:]
pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])

# Custom location, pix size, and image size
# ref_head = make_header(4.08, [256,256], 0.0756, -21.5451)
# pixsize = 4.08

PSW_I_map, ra, dec =  create_map(ref_head, nu=nu.value)
gridra = np.reshape(ra.value,PSW_I_map.shape)
griddec = np.reshape(dec.value,PSW_I_map.shape)

mid_ra = gridra[np.int(len(gridra)/2),np.int(len(gridra)/2)]
mid_dec = griddec[np.int(len(griddec)/2),np.int(len(griddec)/2)]
PSW_header = make_header(pixsize, PSW_I_map.shape, mid_ra, mid_dec)

hdu1 = fits.PrimaryHDU(PSW_I_map, PSW_header)
hdul1 = fits.HDUList([hdu1])
hdul1.writeto('iris_sfd_' + timestamp + '_fx.fits',overwrite=True)
# hdul1.writeto('iris_sfd_test.fits',overwrite=True)

print('done')

# hdu2 = fits.PrimaryHDU(gridra, PSW_header)
# hdul2 = fits.HDUList([hdu2])
# hdul2.writeto('planck_' + timestamp + '_ra.fits',overwrite=True)
#
# hdu3 = fits.PrimaryHDU(griddec, PSW_header)
# hdul3 = fits.HDUList([hdu3])
# hdul3.writeto('planck_' + timestamp + '_dc.fits',overwrite=True)

# min_dec = np.min(dec)
# max_dec = np.max(dec)
# min_ra  = np.min(ra)
# max_ra  = np.max(ra)
#
# fig,ax = plt.subplots(figsize=(14,11))
# plt.imshow(PSW_I_map, origin='lower', extent=[min_dec, max_dec, min_ra, max_ra])#, clim=(1.8916812, 8.812404))
# cbar = plt.colorbar()
# plt.xlabel('Dec')
# plt.ylabel('RA')
# cbar.set_label(r'$I_{\nu}$ [MJy/sr]')
# plt.show()
