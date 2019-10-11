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
from mosaiq.py import mosaic
import numpy as np
from astropy.io import fits
from utility import *
from astropy.wcs import WCS as world

#choose the field that you want
field = config.fields[4]
ra = field[0]
dec = field[1]

#default pixel size of LORRI
pixsize = 4.1 / 3600 #pixel size in arcseconds / pixel

#want a large path of the sky (iris images are 500, 500)
naxis = int(sqrt(2.) * 256)

#square patch of sky, this can be changed if desired.
naxis1 = naxis2 = naxsis

#holder for the x, y pixel coordinates that we want.
xvec = np.arange(0, naxis1)
yvec = np.arange(0, naxis2)

head = make_header(pixsize, naxis, ra, dec) #this function assumes a square.

w = world(head)
c = pixel_to_skycoord(xvec, yvec, w, origin=0)
#this is converting the pixel coords to right ascension and declination in fk4
ra = np.asarray(c.ra.to_string(decimal=True), dtype=float)
dec = np.asarray(c.dec.to_string(decimal=True), dtype=float)

#get information from the IRAS astrometry positions
map = mosaic(header, band=4)
