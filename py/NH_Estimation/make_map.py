################################################################################
# NAME : make_map.py
# DATE STARTED : June 25, 2020
# AUTHORS : Benjamin Vaughan
# PURPOSE : the purpose of this program is to make a map of Cirrus Emission using
# the results of the Dust Fitting Analysis done by the Planck Team on their 2015
# GNILC-Dust models.
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
from astropy.io import fits
import numpy as np
from math_functions import *
from scipy.interpolate import griddata
from astropy.wcs import WCS as world
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import matplotlib.pyplot as plt
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import griddata
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
import nh_config
from reproject import reproject_interp, reproject_from_healpix

def create_map(ref_head, nu):
    '''
    Inputs : filenames (str array) - the paths to and the filenames of the component maps in order of Tau, Temp, Beta (see example.py)
             ref_head  (astropy.header) - the header for the reference map to interpolate onto
             ref_pixsize (float) - the pixel size of the reference map
             ref_mapsize (int array) - a set of the values for the length of each axes of the map
             center (float array) - a set of the values for the center of the field of interest
    Outputs : Interped_map (float array) - an array of fluxes interpolated onto the reference map's grid.
              ragrd (float array) - an array of the Right Ascension values used in the interpolation
              decgrd (float array) - an array of the Declination values used in the interpolation
    Notes : The Planck model has an offset in intensity due to the fact that in Planck Collaborationet al. 2016 correlated their
            fit to regions where H1 is present and therefore had an offset in their fit. For more details see their paper.
    '''

    filenames = ['NHI_HPX.fits']

    #fix astrometry assuming square (astropy prefers CD matrices)
    if 'CD1_1' in ref_head.keys():
        ref_pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])
    elif 'CDELT1' in ref_head.keys():
        ref_head['CD1_1'] = - 1 * ref_head['CDELT1']
        ref_head['CD1_2'] = 0
        ref_head['CD2_1'] = 0
        ref_head['CD2_2'] = ref_head['CDELT2']
        ref_pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])

    # I_to_MJy = 1e20 #converts from standard units of Intensity to MJy
    param_values = []
    for f in filenames:
        data, ra_grid, dec_grid = read_in_fits(f, ref_head)
        param_values.append(data)
    # I = calc_intensity(param_values[2], param_values[0], param_values[1], nu)

    I_map = param_values[0]

    return I_map, ra_grid, dec_grid
    # return interped_map

def read_in_fits(filename, ref_head):
    '''
    Purpose : This function reads in the fits files for the components and parses them so that we are left with data only for our field.
    Inputs: filename (str) - the singular filename for a component used in the MBB fit
            ref_head (Astropy.header) - header for the reference field
    Outputs: map (float array) - an array of flux values at the given field
             RA_grid (float array) - grid of the Right Ascension values used to pull out components
             DEC_grid (float array) - grid of the Declination values used to pull out components
    '''

    hdul = fits.open(filename)
    head = hdul[1].header
    data = hdul[1].data.field(5) #6th column has map NH data
    nside = head['NSIDE']
    order = head['ORDERING']
    hdul.close()

    if order == 'RING':
        nested = False
    else:
        nested = True

    map, mask = reproject_from_healpix(filename, world(ref_head), shape_out=(ref_head['NAXIS1'], ref_head['NAXIS2']), nested=nested, field=5)
    # print(map.shape)
    x  = np.arange(0, ref_head['NAXIS2'])
    y  = np.arange(0, ref_head['NAXIS1'])

    X, Y = np.meshgrid(x, y)
    w = world(ref_head)
    skycoords = pixel_to_skycoord(X.ravel(), Y.ravel(), wcs=world(ref_head), origin=0)
    RA_grid = np.asarray(skycoords.ra.to_string(decimal=True), dtype='float') * u.deg
    DEC_grid = np.asarray(skycoords.dec.to_string(decimal=True), dtype='float') * u.deg
    return map, RA_grid, DEC_grid


def interp_back_to_ref(img, ra, dec, ref_head):
    '''
    this function is depreciated use reproject_interp from astropy instead
    Purpose: Perform the inperolation of the completed Planck Map to a reference header
    Inputs: img - the Planck flux map
            ra  - grid of the Right Ascension values used to create the Planck Map
            dec - grid of the Decliation values used to create the Planck Map
            ref_head - the reference header for the interpolatino grid
            ref_shape - the shape of the interpolation grid
    '''

    ref_shape = [ref_head['naxis2'], ref_head['naxis1']]

    map_size = img.shape

    # reformat map data and coordinates
    data = np.ravel(img)
    points = np.column_stack((np.ravel(ra), np.ravel(dec)))

    #create a agrid of RA/DEC coordinates that we want to interpolate over
    ref_w = world(ref_head)
    ref_grid_x, ref_grid_y = np.mgrid[0:ref_shape[0], 0:ref_shape[1]]
    ref_grid_ra, ref_grid_dec = ref_w.wcs_pix2world(ref_grid_x, ref_grid_y, 0)

    #do the interpolation
    interp_map = griddata(points, data, (ref_grid_ra, ref_grid_dec), method='linear')
    final_map = np.swapaxes(interp_map, 0, 1)
    return final_map, ref_grid_ra, ref_grid_dec
