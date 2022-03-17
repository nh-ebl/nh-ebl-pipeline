
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
import sys
import os
from mosaic import mosaic
import numpy as np
import config
from astropy.io import fits
from utility import *
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import FK4

#choose the field that you want
#fieldnum is the label of the field
#fieldnum = 1
#num corresponds to the field number in the fields.txt list and not related to the matlab label, used to index fields.txt
#field = config.fields[:,fieldnum] # : selects both RA and Dec
currentDir = os.path.dirname(os.path.realpath(__file__)) #get current working directory (calling from cmd line python has directory issues)
os.chdir(currentDir) #change dir to the needed directory
fileID = open("imagefile.txt","r") #open the file
imagefile = fileID.readlines()[0] #get the imagefile path
fileID.close() #close the file
hdul = fits.open(imagefile) #use the imagefile path
#this is a reference image
# hdul = fits.open('regist_pluto_20151102_030873_lor_0308731028_0x633_sci.fit')

ref_head = hdul[0].header
timestamp = ref_head['SPCUTCJD'][3:]
pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])

# # Retrieve image with specified ra/dec, pixsize in degrees, and image width (assuming square)
# pixsize = 4.08
# ref_head = make_header(pixsize/3600, 256, 0.0756, -21.5451)

#note that the fields in fields.txt are given in the ICRS coordinate frame, but the IRAS data has B1950 Astrometry.
# ra_c = ref_head['SPCBRRA'] *u.deg
# dec_c = ref_head['SPCBRDEC'] * u.deg
# coord = SkyCoord(ra=ra_c, dec=dec_c)
# b1950_coord = coord.transform_to(FK4(equinox='B1950'))
# ra_c_b = float(b1950_coord.ra.to_string(decimal=True))
# dec_c_b = float(b1950_coord.dec.to_string(decimal=True))

naxis = 256
#square patch of sky, this can be changed if desired.
naxis1 = naxis2 = naxis

# y = np.arange(0, naxis1)
# x = np.arange(0, naxis2)
# yvec = np.repeat(x[:, np.newaxis], naxis2, axis=1)
# xvec = np.repeat(y[np.newaxis, :], naxis1, axis=0)
# head = make_header(pixsize, naxis, ra_c_b, dec_c_b) #this function assumes a square.

# w = world(head)
# c = pixel_to_skycoord(xvec, yvec, w, origin=0)
# c = c.transform_to('icrs')

#this is converting the pixel coords to right ascension and declination in fk4
# ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
# dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)

#get information from the IRAS astrometry positions
map = mosaic(ref_head, band=4)

hdu = fits.PrimaryHDU(map, ref_head)
hdul = fits.HDUList([hdu])
hdul.writeto(config.DataDir + 'iris_' + timestamp + '_fx.fits',overwrite=True)
# hdul.writeto('iris_test.fits',overwrite=True)

print('done')

#ra
# hdu = fits.PrimaryHDU(ra, head)
# hdul = fits.HDUList([hdu])
# hdul.writeto(config.DataDir + 'iris_' + timestamp + '_ra.fits')
# #dec
# hdu = fits.PrimaryHDU(dec, head)
# hdul = fits.HDUList([hdu])
# hdul.writeto(config.DataDir + 'iris_' + timestamp + '_dc.fits')