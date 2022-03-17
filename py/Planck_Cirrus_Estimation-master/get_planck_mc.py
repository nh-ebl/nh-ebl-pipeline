################################################################################
# NAME : model_test.py
# DATE STARTED : July 5, 2020
# AUTHORS : Benjamin Vaughan, Teresa Symons
# PURPOSE : The purpose of this script is to test that the map we recreate using
# the best fit params from Planck agree with their model.
# Run it N times for N means of the image instead of a single image
# EXPLANATION :
# CALLING SEQUENCE :
# INPUTS :
#
#
# OUTPUTS :
# REVISION HISTORY :
################################################################################
import sys
import os
from astropy.io import fits
sys.path.append('../')
from utilities import *
from make_map import *
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.wcs import WCS as world
from astropy import units as u
from astropy.constants import c
from scipy import io
from joblib import Parallel, delayed
from multiprocessing import cpu_count

# Flag to tell make_map if varying parameters by error or not
FLG_err = 1

# Number of times to run this script, hard-coded
N = 100
seed = 1776

lam = 100*u.micron
nu = (c.to(u.micron/u.s)/lam).to(u.Hz)

#this is a reference image
# hdul = fits.open('regist_pluto_20151102_030873_lor_0308731028_0x633_sci.fit')
#<read in the fits image name from a text file made by matlab>
currentDir = os.path.dirname(os.path.realpath(__file__)) #get current working directory (calling from cmd line python has directory issues)
os.chdir(currentDir) #change dir to the needed directory
fileID = open("imagefile.txt","r") #open the file
imagefile = fileID.readlines()[0] #get the imagefile path
fileID.close() #close the file
hdul = fits.open(imagefile) #use the imagefile path


ref_head = hdul[0].header
timestamp = ref_head['SPCUTCJD'][3:]
pixsize = 3600 *  np.mean([abs(ref_head['CD1_1'] + ref_head['CD2_1']), abs(ref_head['CD2_1'] + ref_head['CD2_2'])])

map_err_mean = {'Temperature':np.zeros(N), 'Spectral-Index':np.zeros(N), 'Opacity':np.zeros(N)}
FLG_noerrPrecalcd = create_map(ref_head, FLG_err, seed, nu=nu.value,FLG_noerrPrecalcd=True) #load in no error I_maps so they do not need to be recalculated
seed_array = seed + np.arange(0,100,1) #precalculate
#--- non-parallel method ---
# for i in range(0,N):
#     print('On MC simulation: '+str(i))
#     PSW_I_map = create_map(ref_head, FLG_err, seed_array[i], nu=nu.value,FLG_noerrPrecalcd=FLG_noerrPrecalcd)
#     keyz = list(PSW_I_map.keys()) #the order of the keys is not guaranteed to be consistent, so can't rely on the order
#     for j in range(0,len(keyz)):
#         map_err_mean[keyz[j]][i] = np.mean(PSW_I_map[keyz[j]])
#     #seed += 1
# map_err_mean_out = np.zeros((3,N))  #Tau = 0, Temperature = 1, Beta = 2 | Opacity = 0, Temperature = 1, Spectral-Index = 2
# map_err_mean_out[0,:] = map_err_mean['Opacity']
# map_err_mean_out[1,:] = map_err_mean['Temperature']
# map_err_mean_out[2,:] = map_err_mean['Spectral-Index']
#--- parallel method ---
def mc_sim_parallel(ref_head, FLG_err, seed_array_value, nu_value, FLG_noerrPrecalcd, i):
    print('On MC simulation: '+str(i))
    map_err_mean = {'Temperature':0., 'Spectral-Index':0., 'Opacity':0.}
    PSW_I_map = create_map(ref_head, FLG_err, seed_array_value, nu=nu_value,FLG_noerrPrecalcd=FLG_noerrPrecalcd)
    keyz = list(PSW_I_map.keys()) #the order of the keys is not guaranteed to be consistent, so can't rely on the order
    for j in range(0,len(keyz)):
        map_err_mean[keyz[j]] = np.mean(PSW_I_map[keyz[j]])
    #seed += 1
    return map_err_mean

mc_sim_parallel_list = [None for i in range(0,N)] #preallocate list
for i in range(0,N): #create input list for parallel stuff
    mc_sim_parallel_list[i] = [ref_head, FLG_err, seed_array[i], nu.value, FLG_noerrPrecalcd, i]

mc_sim_parallel_results = Parallel(n_jobs=cpu_count())(delayed(mc_sim_parallel)(i, j, k, l, m, n) for i, j, k, l, m, n in mc_sim_parallel_list)

map_err_mean_out = np.zeros((3,N))  #Tau = 0, Temperature = 1, Beta = 2 | Opacity = 0, Temperature = 1, Spectral-Index = 2
for i in range(0,N):
    map_err_mean_out[0,i] = mc_sim_parallel_results[i]['Opacity']
    map_err_mean_out[1,i] = mc_sim_parallel_results[i]['Temperature']
    map_err_mean_out[2,i] = mc_sim_parallel_results[i]['Spectral-Index']

io.savemat('planck_' + timestamp + '_errmean.mat',{'err_means': map_err_mean_out})
print('done')
# gridra = ra
# griddec = dec

# ra  = ra[:,0]
# dec = dec[0, :]
# mid_ra = np.median(ra)
# mid_dec = np.median(dec)
# PSW_header = make_header(pixsize, PSW_I_map.shape, mid_ra, mid_dec)

# hdu1 = fits.PrimaryHDU(PSW_I_map, PSW_header)
# hdul1 = fits.HDUList([hdu1])
#hdul1.writeto('planck_' + timestamp + '_fx.fits',overwrite=True)

# hdu2 = fits.PrimaryHDU(gridra, PSW_header)
# hdul2 = fits.HDUList([hdu2])
# hdul2.writeto('planck_' + timestamp + '_ra.fits',overwrite=True)
#
# hdu3 = fits.PrimaryHDU(griddec, PSW_header)
# hdul3 = fits.HDUList([hdu3])
# hdul3.writeto('planck_' + timestamp + '_dc.fits',overwrite=True)

# min_dec = np.min(dec)
# max_dec = np.max(dec)
# min_ra  = np.min(ra)
# max_ra  = np.max(ra)
#
# fig,ax = plt.subplots(figsize=(14,11))
# plt.imshow(PSW_I_map, origin='lower', extent=[min_dec, max_dec, min_ra, max_ra])#, clim=(1.8916812, 8.812404))
# cbar = plt.colorbar()
# plt.xlabel('Dec')
# plt.ylabel('RA')
# cbar.set_label(r'$I_{\nu}$ [MJy/sr]')
# plt.show()
