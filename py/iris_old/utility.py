################################################################################
# Name : mosaiq
# Purpose : Utility functions that are used when creating an IRIS mosaic
#Author : Benjamin Vaughan
#Start Date : October 4, 2019
#Additional Info
#
################################################################################
from astropy.io import fits
import numpy as np
from scipy.interpolate import griddata
from math import *

def make_header(pixsize, naxis, ra, dec):
    '''
    Purpose : This function makes a fits header with all the astrometry you need to
    do transformations from ra/dec -> pixel and pixel -> ra/dec
    Inputs  : pixsize - (float) the pixel size for your map
              naxis   - (int) the size of your map (this function assumes a square map)
              ra      - (float) right ascension in degrees of the origin
              dec     - (float) declination in degrees of the origin
    Outputs : Header - a fits header object containing astrometry data
    '''
    cd1_1 = -pixsize
    cd1_2 = 0
    cd2_1 = 0
    cd2_2 = pixsize
    crpix = naxis / 2 + 1
    crval1 = ra
    crval2 = dec
    header = fits.Header()
    header.set('NAXIS', int(2))
    header.set('NAXIS1', naxis)
    header.set('NAXIS2', naxis)
    header.set('CD1_1', cd1_1)
    header.set('CD1_2', cd1_2)
    header.set('CD2_1', cd2_1)
    header.set('CD2_2', cd2_2)
    header.set('CRPIX1', crpix)
    header.set('CRPIX2', crpix)
    header.set('CRVAL1', crval1)
    header.set('CRVAL2', crval2)
    header.set('EQUINOX', 'unknown')
    header.set('BITPIX', -64)
    header.set('CTYPE1', 'RA---TAN')
    header.set('CTYPE2', 'DEC--TAN')
    header.set('EQUINOX', 2000)
    return header

def nan2undef(data, undef=-32768):
    '''
    Purpose : Sets bad data to a value of -32768 by default or a specified value
    Inputs  : data  - (2D array) the data you want to check through
              undef - (float) default at -32768, set undefined data to this value
    Outputs : data  - (2D array) data with undefined values set to a specific value
    '''
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if np.isfinite(data[i,j]) != True:
                data[i,j] = undef
    return data

def get_cord_type(header):
    '''
    Purpose : gets the coordinates for an image and tells you what coordinate system
    Inputs : Header - header info
    Outputs: ctype - the coordinate type
    '''
    c1 = header['ctype1'] #functional "x" axis
    c2 = header['ctype2'] #functional "y" axis

    if 'RA' in c1 and 'DEC' in c2:
        ctype = 1 #sky coordinates
    elif 'GLON' in c1 and 'GLAT' in c2:
        ctype = 2 #galactic coordinates
    elif 'ELON' in c1 and 'ELAT' in c3:
        ctype = 3 #Ecliptic              head - the header for the corresponding map

    else:
        print('unknown coordinate system aborint!')
        exit()
    return ctype


def get_array_value(x, y, array):
    """Returns the value of the array at position x,y."""
    return array[y, x]

def linterp(x, y, ratio):
    '''
    The purpose of this function is to honestly make the code look neat. It does
    a linear interpolation across one axis.
    Inputs : x - the x position
             y - the y position
             ratio - an interpolation ratio calculated from the interpolant
                     coordinates
    Returns: interpolant
    '''

    return x * (1 - ratio) + y * ratio

def bilinear_interpolation(x, y, img):
    '''
    This function is meant to perform a bilinear interpolation of an image
    Inputs : x - an array of where you want array to be interpolated at
             y - an array of where you want array to be interpolated at
             array - array to be interpolated
    Returns: result - interpolant
    '''

    result = []
    for i in range(len(x)):
        x2 = int(ceil(x[i]))
        x1 = int(floor(x[i]))
        y2 = int(ceil(y[i]))
        y1 = int(floor(y[i]))

        ratio_x = x[i] - x1
        ratio_y = y[i] - y1

        interpolate_x1 = linterp(img[y1, x1],
                                     img[y1, x2],
                                     ratio_x)
        interpolate_x2 = linterp(img[y2, x1],
                                     img[y2, x2],
                                     ratio_x)
        interpolate_y = linterp(interpolate_x1, interpolate_x2, ratio_y)
        result.append(interpolate_y)

    result = np.asarray(result)

    return result
