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

def interpolate(first_value, second_value, ratio):
    """Interpolate with a linear weighted sum."""
    return first_value * (1 - ratio) + second_value * ratio

def bilinear_interpolation(x, y, img):
    """Returns the bilinear interpolation of a pixel in the image.
    :param x: x-position to interpolate
    :param y: y-position to interpolate
    :param img: image, where the pixel should be interpolated
    :returns: value of the interpolated pixel
    """
    # if x < 0 or y < 0:
    #     raise ValueError('x and y pixel position have to be positive!')
    # if img.shape[1]-1 < x or img.shape[0]-1 < y:
    #     raise ValueError('x and y pixel position have to be smaller than image'
    #                      'dimensions.')

    listy = []
    for i in range(len(x)):
        x_rounded_up = int(ceil(x[i]))
        x_rounded_down = int(floor(x[i]))
        y_rounded_up = int(ceil(y[i]))
        y_rounded_down = int(floor(y[i]))

        ratio_x = x[i] - x_rounded_down
        ratio_y = y[i] - y_rounded_down

        interpolate_x1 = interpolate(img[y_rounded_down, x_rounded_down],
                                     img[y_rounded_down, x_rounded_up],
                                     ratio_x)
        interpolate_x2 = interpolate(img[y_rounded_up, x_rounded_down],
                                     img[y_rounded_up, x_rounded_up],
                                     ratio_x)
        interpolate_y = interpolate(interpolate_x1, interpolate_x2, ratio_y)
        listy.append(interpolate_y)

    listy = np.asarray(listy)

    return listy
