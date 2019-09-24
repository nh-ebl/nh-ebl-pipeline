################################################################################
# Name : mosaiq
# Purpose : Copy of idl MOSAIQ script from the french people in python
# I guess this creates a mosaic of IRIS images?
#Author : Benjamin Vaughan
#Additional Info
#
################################################################################
from astropy.io import fits
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as wcs

def mosaic(header, band=None, catname=None, dir=None):
    '''
    purpose: create mosaic of images
    Inputs : Header - header information
             band   - band (defaults to 4 in idl)
             catname- filename of the catalog (defaults to !IRISPRO which is some config file)
             dir    - directory where iris maps are (default is !IRISDATA which is some config file)
    Outputs:
    '''
    #w = wcs(header) # there are issues with this for some reason.
    try:
        equinox = header['EQUINOX'] #get equinox variable from a fits header #doesn't seem to be an equinox in the fits files?????
    except KeyError:
        equinox = header['EPOCH'] #this is what they did in NASA's idl extast function

    x_size = header['NAXIS1']
    y_size = header['NAXIS2']
    xmap = np.arange(0, x_size)
    ymap = np.arange(0, y_size)

    #ra, dec = pixel_to_skycoord(x_size, y_size, w, origin=0) #origin is either 1 or 0 can't remember which corresponds to what
    ctype = get_cord_type(header)
    if ctype == 1:
        fk4 = 1 #what does fk4 mean aside from the obvious coordinate.

    if ctype == 2:
        fk4 = 1
        equinox = 1950 #force the equinox to be 1950 in the case of galactic coordinates
        print('converting from Galactic to Celestial')
        # do the thing
    elif ctype == 3:
        print('converting from Ecliptic to Celestial')
        # do the thing

    ra = nan2undef(ra)
    dec = nan2undef(dec)

    #open up the info_issa_map4.txt file and read it into things that we can use
    #need to figure out how it reads file if it does anything special, doesn't seem like it...

    #nb = n_elements(inum) this means find the number of columns

    print('Checking for ISSA maps that intersect with the given header')

    ind1 = np.where(ra != -32768)
    ind2 = np.where(dec != -32768)

    combined_ind = []

    good_inds = np.zeros(nb)

    if ind1 in ind2:
        combined_ind.append(ind1)

    c1min = min(ra[combined_ind])
    c1max = max(ra[combined_ind])
    c2min = min(dec[combined_ind])
    c2max = max(dec[combined_ind])

    for i in range(nbind):
        pass #this requires us to do the text file read in

    good_inds = np.where(good_inds > 0)
    if len(good_inds) > 0:
        print('No ISSA map corresponds to the header given')
        exit()

    print('%s ISSA maps will be combined to produce the mosaic' %(nbind))

    for i in range(nbind):
        #do the actual point of this script
        #call get_iris
        pass

    indw = np.where(weight > 0)
    if len(indw) > 0:
        result[indw] = result[indw] / weight[indw]

    badindwmap = np.arange(0, weight.shape)
    badind = []
    for i in range(len(badindw)):
        if indw in badindw:
            pass
        else:
            badind.append(i)

def get_cord_type(header):
    '''
    Purpose : gets the coordinates for an image and tells you what coordinate system
    Inputs : Header - header info
    Outputs: ctype - the coordinate type
    '''
    c1 = header['ctype1']) #functional "x" axis
    c2 = header['ctype2']) #functional "y" axis

    if 'RA' in c1 and 'DEC' in c2:
        ctype = 1 #sky coordinates
    elif 'GLON' in c1 and 'GLAT' in c2:
        ctype = 2 #galactic coordinates
    elif 'ELON' in c1 and 'ELAT' in c3:
        ctype = 3 #Ecliptic
    else:
        print('unknown coordinate system aborint!')
        exit()
    return ctype

def bprecess():
    # needs to be written
    # some preliminary documentation on what to look up for this:
    # https://docs.astropy.org/en/stable/coordinates/index.html#module-astropy.coordinates
    # https://pythonhosted.org/Astropysics/index.html
    # https://idlastro.gsfc.nasa.gov/ftp/pro/astro/jprecess.pro
    # https://www.astrobetter.com/wiki/Python+Switchers+Guide
    print(' precessing coordinates from J2000 to B1950')
    pass #docmentation poitns to jprecess we actually want bprecess

def nan2undef(data, undef=-32768, indef=False):
    if indef:
        undef = indef
    ind = np.where( data != 1)
    if ind > 0:
        data[ind] = undef
    return data





if __name__ == '__main__':
    file = '../IRISNOHOLES_B4H0//I005B4H0.FIT'
    f = fits.open(file)
    for d in f[0].header:
        print(d)
    mosaic(f[0].header)
