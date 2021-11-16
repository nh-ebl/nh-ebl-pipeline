################################################################################
# Name : nh_iris_map
# Purpose : Copy of the IDL code that Mike Zemcov wrote to talk to the Mosaique
# scripts, but ported over to python. It continues multiple potential fields
# depending on what you want the location of the origin for your new header
# information to be. This creates a mosaic with all of the IRIS maps that
# intersect with the given field
#Author : Benjamin Vaughan
#Start Date : Oct 11, 2019
#Additional Info
#
################################################################################
from mosaiq import mosaic
import numpy as np
import config
from astropy.io import fits
from utility import *
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt

#choose the field that you want
#fieldnum is the label of the field, corresponding to how fields are labeled in the matlab code, only used to save file name
fieldnum = 8
#num corresponds to the field number in the fields.txt list and not related to the matlab label, used to index fields.txt
field = config.fields[:,4] # : selects both RA and Dec

ra = field[0]
dec = field[1]

#default pixel size of LORRI
pixsize = 4.1 / 3600 #pixel size in arcseconds / pixel

#want a large path of the sky (iris images are 500, 500)
naxis = int(sqrt(2.) * 256)

#square patch of sky, this can be changed if desired.
naxis1 = naxis2 = naxis

#holder for the x, y pixel coordinates that we want.
y = np.arange(0, naxis1)
x = np.arange(0, naxis2)
yvec = np.repeat(x[:, np.newaxis], naxis2, axis=1)
xvec = np.repeat(y[np.newaxis, :], naxis1, axis=0)
head = make_header(pixsize, naxis, ra, dec) #this function assumes a square.

w = world(head)
c = pixel_to_skycoord(xvec, yvec, w, origin=0)
#this is converting the pixel coords to right ascension and declination in fk4
ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)

#get information from the IRAS astrometry positions
map = mosaic(head, band=4)

#for debugging
# plt.imshow(map, origin='lower')
# plt.show()

hdu = fits.PrimaryHDU(map, head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_0' + str(fieldnum) + '_fx.fits')

#ra
hdu = fits.PrimaryHDU(ra, head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_0' + str(fieldnum) + '_ra.fits')

#dec
hdu = fits.PrimaryHDU(dec, head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_0' + str(fieldnum) + '_dc.fits')
