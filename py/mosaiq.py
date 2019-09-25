################################################################################
# Name : mosaiq
# Purpose : Copy of idl MOSAIQ script from the french people in python
# I guess this creates a mosaic of IRIS images?
#Author : Benjamin Vaughan
#Start Date : September 23, 2019
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

    #there's a better way to do this than defining like 8 variables and setting it equal to this function output but oh well
    inum, ramin, ramax, raavg, decmin, decmax, decavg, medianval, noise_key = np.loadtxt(catname, unpack=True)
    #nb = n_elements(inum) this means find the number of columns
    numel = inum[-1] #number of entries in the text file
    print(numel)
    print('Checking for ISSA maps that intersect with the given header')

    ind1 = np.where(ra != -32768)
    ind2 = np.where(dec != -32768)

    combined_ind = []

    good_inds = np.zeros(numel)

    if ind1 in ind2:
        combined_ind.append(ind1)

    c1min = min(ra[combined_ind])
    c1max = max(ra[combined_ind])
    c2min = min(dec[combined_ind])
    c2max = max(dec[combined_ind])


    #i feel as if a lot of these checks are redundant, but there is no need to go to deep into figuring out
    #a better way to do this, because we don't really gain anything from that unless we are going to run this
    #literally a million times
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

    good_inds = np.where(good_inds > 0)
    if len(good_inds) > 0:
        print('No ISSA map corresponds to the header given')
        exit()

    print('%s ISSA maps will be combined to produce the mosaic' %(nbind))

    for i in range(nbind):
        mapi = get_iris(inum[i], dir=dir, band=band)
        #converting alpha and delta back to pixel coords?
        #if all alpha and delta is doing before this is finding
        #min / max RAs and Decs there is no need to convert the entire
        #array to RA and Dec because that is stupid and would take forever
        tempo = mbillinear() #this function needs to be written 

    indw = np.where(weight > 0)
    if len(indw) > 0:
        result[indw] = result[indw] / weight[indw]

    badindwmap = np.arange(0, weight.shape)
    badind = []
    for i in range(len(badindw)):
        if indw in badindw:

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

def get_iris(ii, dir='!ISRISDATA', band=4, hcon=0, verbose=0):
    '''
    Purpose : returns the correlated IRIS map for an input ISSA number
    Inputs  : ii - the ISSA map number
              verbose - verbosity flag
              band - The IRAS band
                - 1 : 12 micron
                - 2 : 25 micron
                - 3 : 60 micron
                - 4 : 100 micron (default)
              hcon - hcon number (default is zero -> co-added map)
              dir - directory where the IRIS data is stored (default is !IRISDATA in idl)
    Outputs : map - the ISSA map corresponding to number ii
    '''
    iras_number = str(ii)
    if ii < 10:
        iras_number = '0' + iras_number
    elif ii < 100:
        iras_number = '0' + iras_number # um maybe there's a difference in idl but this if statement makes no sense

    files = []
    for x in os.listdir(dir):
        if iras_number in x and bd in x and hcon in x:
            file.append(x) #this needs to be tested more rigourously, but im sure this is right
    if len(files) > 0:
        hdul = fits.open(file[0])
        header = hdul[0].header
        header.append(('LONPOLE',180))
        bad = np.where(hdul[0].data <= -5 or hdul[0].data == 0) #have to verify this method for 2 dimensional arrays
        hdul[0].data[bad] = -32768
        hdu = fits.PrimaryHDU(hdul[0].data, hdul[0].header)
        map = fits.HDUList([hdu])
        if verbose:
            print('Read data file %s' file[0])
        return map
    else:
        print('Could not find any files matching that description')
        return None

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
