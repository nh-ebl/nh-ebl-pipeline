################################################################################
# Name : mosaiq
# Purpose : copy of the IDL routine mosaique_iris.pro ported over to python
# this code does a transformation on an image to rotate and align it with
# another image.
#Author : Benjamin Vaughan
#Start Date : September 23, 2019
#Additional Info
#
################################################################################
from astropy.io import fits
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK4
from astropy.coordinates import Galactic
from astropy.coordinates import BarycentricTrueEcliptic
import astropy.units as u
import sys
from math import *
import os
from bilinear import mbilinear
from utility import *
from get_iris_map import get_iris
import matplotlib.pyplot as plt

def mosaic(header, band=4, catname=None, dir=None):
    '''
    purpose: create mosaic of images
    Inputs : Header - header information
             band   - band (defaults to 4 in idl)
             catname- filename of the catalog (defaults to !IRISPRO which is some config file)
             dir    - directory where iris maps are (default is !IRISDATA which is some config file)
    Outputs:
    '''


    header.set('NAXIS', 2)
    w = world(header)
    try:
        equinox = header['EQUINOX'] #get equinox variable from a fits header #doesn't seem to be an equinox in the fits files?????
    except KeyError:
        equinox = header['EPOCH'] #this is what they did in NASA's idl extast function


    x_size = header['NAXIS1']
    y_size = header['NAXIS2']
    print('Handling %s elements' % (x_size * y_size))
    x = np.arange(0, x_size)
    y = np.arange(0, y_size)
    xlist = np.tile(1, x_size)
    ylist = np.tile(1, y_size)
    xmap = np.outer(x, ylist)
    ymap = np.outer(xlist, y)

    result = np.zeros((x_size, y_size))
    weight = result
    new_c = pixel_to_skycoord(xmap, ymap, w, origin=0)
    #this is converting the pixel coords to right ascension and declination in fk4
    ra = np.asarray(new_c.ra.to_string(decimal=True), dtype=float)
    dec = np.asarray(new_c.dec.to_string(decimal=True), dtype=float)

    ctype = get_cord_type(header)

    if ctype == 2:
        fk4 = 1
        equinox = 1950 #force the equinox to be 1950 in the case of galactic coordinates
        print('converting from Galactic to Celestial')
        sk = SkyCoord(ra * u.deg, dec * u.deg, frame=Galactic)
        new_c = sk.transform_to(FK4)
        ra = np.asarray(new_c.ra.to_string(decimal=True), dtype=float)
        dec = np.asarray(new_c.dec.to_string(decimal=True), dtype=float)

    elif ctype == 3:
        print('converting from Ecliptic to Celestial')
        #need to know what ecliptic it is.
        sk = SkyCoord(ra * u.deg, dec * u.deg, frame=BarycentricTrueEcliptic)
        new_c = sk.transform_to(FK4)

        ra = []
        dec = []
        coordinates = new_c.to_string('decimal')
        ra = np.asarray(new_c.ra.to_string(decimal=True), dtype=float)
        dec = np.asarray(new_c.dec.to_string(decimal=True), dtype=float)

    if equinox == 2000.0:
        print('precessing coordinates from J2000 to B1950')
        new_c.transform_to(FK4(equinox='B1950'))
        ra = np.asarray(new_c.ra.to_string(decimal=True), dtype=float)
        dec = np.asarray(new_c.dec.to_string(decimal=True), dtype=float)

    ra = nan2undef(ra)
    dec = nan2undef(dec)
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    #converted these to arrays for easier data manipulation

    inum, ramin, ramax, raavg, decmin, decmax, decavg, medianval, noise_key = np.loadtxt(catname, unpack=True)
    numel = int(inum[-1]) #number of entries in the text file
    print('Checking for ISSA maps that intersect with the given header')

    ind1 = np.where(ra != -32768)[0]
    ind2 = np.where(dec != -32768)[0]

    ind = []
    for i in range(ra.shape[0]):
        for j in range(ra.shape[1]):
            if ra[i,j] != -32768 and dec[i,j] != -32768:
                ind.append([i,j])
    ind = np.asarray(ind)
    good_inds = np.zeros(numel)

    c1min = np.min(ra[ind])
    c1max = np.max(ra[ind])
    c2min = np.min(dec[ind])
    c2max = np.max(dec[ind])

    for i in range(numel):
        if c1min > ramin[i] and c1min < ramax[i] and c2min > decmin[i] and c2min < decmax[i]:
            good_inds[i] = 1
        elif c1min > ramin[i] and c1min < ramax[i] and c2max > decmin[i] and c2max < decmax[i]:
            good_inds[i] = 1
        elif c1max > ramin[i] and c1max < ramax[i] and c2max > decmin[i] and c2max < decmax[i]:
            good_inds[i] = 1
        elif c1max > ramin[i] and c1max < ramax[i] and c2min > decmin[i] and c2min < decmax[i]:
            good_inds[i] = 1
        elif ramin[i] > c1min and ramin[i] < c1max and decmin[i] > c2min and decmin[i] < c2min:
            good_inds[i] = 1
        elif ramax[i] > c1min and ramax[i] < c1max and decmin[i] > c2min and decmin[i] < c2min:
            good_inds[i] = 1
        elif ramin[i] > c1min and ramin[i] < c1max and decmax[i] > c2min and decmax[i] < c2min:
            good_inds[i] = 1
        elif ramax[i] > c1min and ramax[i] < c1max and decmax[i] > c2min and decmax[i] < c2min:
            good_inds[i] = 1

    #now that we have the maps we want do the outer product trick and give us a 2D array for x and y


    good_inds = np.where(good_inds > 0)[0]
    if good_inds.shape[0] == 0:
        print('No ISSA map corresponds to the header given')
        exit()

    print('%s ISSA maps will be combined to produce the mosaic' %(good_inds.shape[0]))

    for i in range(good_inds.shape[0]):
        if i == 2:
            mapi = get_iris(inum[good_inds[i]], dir=dir, band=band)

            #forcing it to be 2D rather than 3D
            mapi[0].header['NAXIS'] = 2
            mapi[0].header['EPOCH'] = 2000.0
            try:
                del(mapi[0].header['NAXIS3'])
            except KeyError:
                pass

            #do the transform back to pixel coords
            w = world(mapi[0].header)
            x, y = skycoord_to_pixel(new_c, w, origin=0)
            tempo = mbilinear(x, y, mapi[0].data)
            indw = []
            for j in range(tempo.shape[0]):
                for k in range(tempo.shape[1]):
                    if tempo[j,k] != -32768:
                        indw.append([j,k])
            indw = np.asarray(indw)
            weight[indw] = weight[indw] + 1
            result[indw] = result[indw] + tempo[indw]
    indw = []
    complement = []
    for i in range(weight.shape[0]):
        for j in range(weight.shape[1]):
            if weight[i,j] > 0:
                indw.append([i,j])
            else:
                complement.append([i,j])
    complement = np.asarray(complement)
    indw = np.asarray(indw)
    if len(indw) > 0:
        result[indw] = result[indw] / weight[indw]

    result[complement] = -32768

    plt.imshow(result, origin='lower')
    plt.show()
    return result



if __name__ == '__main__':
    #test code to generate input header from Mike's Script
    f = fits.open('../../../../IRISNOHOLES_B4H0/I088B4H0.FIT')
    pixsize = 4.1 / 3600.
    naxis = int(sqrt(2.) * 256)
    # naxis = 500
    naxis1 = naxis2 = naxis
    #the repvec function just replicates vectors
    xvec = np.arange(0, naxis1)
    yvec = np.arange(0, naxis2)
    ra = 196.0075
    dec = 23.9506
    header = make_header(pixsize, naxis, ra, dec)
    catfile = 'info_issa_map4.txt'
    dir = '../../../../IRISNOHOLES_B4H0'
    mosaic(header, catname=catfile, dir=dir)
