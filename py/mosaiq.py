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
from astropy.wcs import WCS as world
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK4
from astropy.coordinates import Galactic
from astropy.coordinates import BarycentricTrueEcliptic
import astropy.units as u
import sys
from math import *
import os
# from mpl_toolkits.basemap import interp
np.set_printoptions(threshold=sys.maxsize)
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
    xmap = np.arange(0, x_size)
    ymap = np.arange(0, y_size)
    xlist = np.tile(1, x_size)
    ylist = np.tile(1, y_size)


    result = np.zeros((x_size, y_size))
    weight = result
    new_c = pixel_to_skycoord(xmap, ymap, w, origin=1) #origin is either 1 or 0 can't remember which corresponds to what
    #this is converting the pixel coords to right ascension and declination in fk4
    ra = []
    dec = []
    coordinates = new_c.to_string('decimal')
    for i in range(len(coordinates)):
        split_coords = coordinates[i].split(' ')
        ra.append(float(split_coords[0]))
        dec.append(float(split_coords[1]))

    ctype = get_cord_type(header)

    if ctype == 2:
        fk4 = 1
        equinox = 1950 #force the equinox to be 1950 in the case of galactic coordinates
        print('converting from Galactic to Celestial')
        sk = SkyCoord(ra * u.deg, dec * u.deg, frame=Galactic)
        new_c = sk.transform_to(FK4)
        ra = []
        dec = []
        coordinates = new_c.to_string('decimal')

        for i in range(len(coordinates)):
            split_coords = coordinates[i].split(' ')
            ra.append(float(split_coords[0]))
            dec.append(float(split_coords[1]))

    elif ctype == 3:
        print('converting from Ecliptic to Celestial')
        #need to know what ecliptic it is.
        sk = SkyCoord(ra * u.deg, dec * u.deg, frame=BarycentricTrueEcliptic)
        new_c = sk.transform_to(FK4)

        ra = []
        dec = []
        coordinates = new_c.to_string('decimal')

        for i in range(len(coordinates)):
            split_coords = coordinates[i].split(' ')
            ra.append(float(split_coords[0]))
            dec.append(float(split_coords[1]))

    if equinox == 2000.0:
        print('precessing coordinates from J2000 to B1950')
        new_c.transform_to(FK4(equinox='B1950'))
        ra = []
        dec = []
        coordinates = new_c.to_string('decimal')

        for i in range(len(coordinates)):
            split_coords = coordinates[i].split(' ')
            ra.append(float(split_coords[0]))
            dec.append(float(split_coords[1]))

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

    combined_ind = []
    for i in range(len(ind1)):
        if ind1[i] in ind2:
            combined_ind.append(ind1[i])

    combined_ind = np.asarray(combined_ind)
    combined_ind = np.unique(combined_ind)

    good_inds = np.zeros(numel)

    c1min = np.min(ra[combined_ind])
    c1max = np.max(ra[combined_ind])
    c2min = np.min(dec[combined_ind])
    c2max = np.max(dec[combined_ind])

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
        if i == 5:
            mapi = get_iris(inum[good_inds[i]], dir=dir, band=band)
            w = world(mapi[0].header)
            xmap, ymap = skycoord_to_pixel(new_c, w, origin=1)
            print(xmap, ymap)
            x = np.outer(xmap, ylist)
            y = np.outer(xlist, ymap)
            tempo = mbilinear(x, y, mapi[0].data) #this function needs to be written
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
    return result

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
    iras_number = str(int(ii))
    if ii < 10:
        iras_number = '0' + iras_number
    elif ii < 100:
        iras_number = '0' + iras_number # um maybe there's a difference in idl but this if statement makes no sense

    files = []
    for x in os.listdir(dir):
        if iras_number in x and str(band) in x and str(hcon) in x:
            files.append(dir + '/' + x)
    if len(files) > 0:
        hdul = fits.open(files[0])

        #fix the header up a little
        header = hdul[0].header
        header.append(('LONPOLE',180))
        if header['CDELT3'] == 0 and header['NAXIS3'] == 1:
            del(header['CDELT3']) #having a value of zero for CDELT3 messes with the convserion
            del(header['NAXIS3']) #if the sizde of the 3rd axis is 1 what is the point of including it.
            # it just messes up the coordinate conversions
            del(header['CRPIX3'])
            del(header['CRVAL3'])
            del(header['CTYPE3'])
            #all of these values will mess up our conversion at the end
            header.set('NAXIS', 2)
            print('Depreciated CDELT3 and NAXIS3 values found, removing.')

        #check data for bad values
        for i in range(hdul[0].data.shape[0]):
            for j in range(hdul[0].data.shape[1]):
                #we have to add this extra zero because for some reason the IRIS maps are
                #3 Dimensional with a z axis that is 1 element thick.
                if hdul[0].data[i,j,0] < -5 or hdul[0].data[i,j,0] == 0:
                    hdul[0].data[i,j] = -32768
        hdu = fits.PrimaryHDU(hdul[0].data, hdul[0].header)
        map = fits.HDUList([hdu])
        if verbose:
            print('Read data file %s' % file[0])
        return map
    else:
        print('Could not find any files matching that description')
        return None

def mbilinear(x, y, array):
    six = x.shape
    siy = y.shape
    sia = array.shape

    #this function has some weird things and it says x and y should be 2D arrays
    Nax = sia[0]
    Nay = sia[1]
    Nx = six[0]
    Ny = six[1] #confused on this because in the idl script it says Ny is the second dim of x...
    #i.e. it would be six[1] but this doesn't make sense to me and the way x and y are generated
    #in the idl version is that they are 1D arrays...
    output = np.zeros((Nx, Ny)) #here it says & output(*,*) = missing, I don't know a python equivalent for this
    #or even what this is really doing (setting every value to missing?)
    #then it finds indexees where the array isn't missing data I guess this would be
    #where the array isn't nan?
    minval = np.min(array)
    # if len(ind) > 0: #i think this is just checking for missing data
    #     minval = np.min(array[ind])

    count = 0
    for i in range(Nx):
        for j in range(Ny):
            if x[i,j] < 0 or x[i,j] > Nax-1 or y[i,j] < 0 or y[i,j] > Nay-1:
                count += 1
    inter_percent = 1. * (Nx * Ny - count) / (Nx * Ny) * 100
    print('Images Intersection = %s percent' % (inter_percent))

    #this line makes no sense when x and y are 1D arrays, I guess this is where the z axis is important :(
    for i in range(Ny):
        mind = []
        for j in range(Nx):
            if x[j,i] >= 0 and x[j,i] <= Nax-1 and y[j,i] >= 0 and y[j,i] <= Nay-1:
                mind.append(j)
        mind = np.asarray(mind)
        if len(mind) != 0:
            xx = np.zeros((len(mind), 2))
            xx[:,0] = x[mind,i]
            yy = np.zeros((len(mind), 2))
            yy[:,0] = y[mind,i]
            # truc = interp(array, xin=xx[:,0], yin=yy[:,0]) #need an interpolation function
            output[ind, i] = truc[:,0] #maybe this needs to be changed.

    #remove values affected by indef values (for highly < 0 indef and generaly > 0 im)
    ind = []
    for i in range(output.shape[0]):
        for j in range(output.shape[1]):
            ind.append([i,j])
    ind = np.asarray(ind)
    output[ind] = np.nan #i guess nan is the same as missing but not sure
    return output



def nan2undef(data, undef=-32768, indef=False):
    if indef:
        undef = indef
    ind = np.where( data != 1)[0]
    if ind > 0:
        data[ind] = undef
    return data

def make_header(pixsize, naxis, ra, dec):
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

if __name__ == '__main__':
    #test code to generate input header from Mike's Script
    pixsize = 4.1 / 3600.
    naxis = int(sqrt(2.) * 256)
    naxis1 = naxis2 = naxis
    #the repvec function just replicates vectors
    xvec = np.arange(0, naxis1)
    yvec = np.arange(0, naxis2)
    ra = 196.0075
    dec = 23.9506
    header = make_header(pixsize, naxis, ra, dec)

    # #now we convert some arrays to ra/dec
    # w = wcs(header)
    # c = pixel_to_skycoord(x, y, w, origin=0)
    # ra = []
    # dec = []
    # coordinates = c.to_string('decimal')
    # for i in range(len(coordinates)):
    #     split_coords = coordinates[i].split(' ')
    #     ra.append(float(split_coords[0]))
    #     dec.append(float(split_coords[1]))

    catfile = 'info_issa_map4.txt'
    dir = '../../../IRISNOHOLES_B4H0'
    mosaic(header, catname=catfile, dir=dir)
