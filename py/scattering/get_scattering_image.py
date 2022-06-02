import numpy as np
import matplotlib.pyplot as plt
from astropy_healpix import HEALPix
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic, ICRS
from scipy import interpolate
from scipy.io import savemat
from pandas import read_csv # Faster than numpy's loadtxt by ~max_ang sec
import sys, os
from multiprocessing import Pool
sys.path.append('../')

print('Loading FITS')
#this is a reference image
# hdul = fits.open('regist_pluto_20160714_033078_lor_0330787618_0x633_sci.fit')
#<read in the fits image name from a text file made by matlab>
currentDir = os.path.dirname(os.path.realpath(__file__)) #get current working directory (calling from cmd line python has directory issues)
os.chdir(currentDir) #change dir to the needed directory
fileID = open("imagefile.txt","r") #open the file
imagefile = fileID.readlines()[0] #get the imagefile path
fileID.close() #close the file
# Read in saved conversion factor
fileID = open("/home/symons/nh_ebl_pipeline/matlab/conv.txt","r") #open the file
sbconv = np.float64(fileID.readlines()[0])
fileID.close() #close the file
hdul = fits.open(imagefile) #use the imagefile path

# Read in header info
ref_head = hdul[0].header
timestamp = ref_head['SPCUTCJD'][3:]
ra_cent = ref_head['SPCBRRA']
dec_cent = ref_head['SPCBRDEC']
cent = SkyCoord(ra=ra_cent*u.degree,dec=dec_cent*u.degree,frame='icrs')

print('Loading PSF')
# Read in lauer psf fig data
lauer_psf = read_csv('lauer_psf.csv', delimiter=',', header=None).to_numpy()
# Save angular distance from center of FOV in degrees
ang_dist = lauer_psf[:,0]
# Save DN/s/pix for V=0 star
brightness = lauer_psf[:,1]

print('Loading map')
# Read in Masana new radiance map info - last column has extra column with '5' after some entries
rad_map_new = read_csv('RadianceOut_new.csv', delimiter=',').to_numpy()
# Save G-band luminosity in W m^-2 sr^-1
L_G_new = rad_map_new[:,10]
L_G_new = L_G_new.astype('float64')

# This plot shows that the nested healpix scheme is being used
# plt.scatter(l_gal_new,b_gal_new,s=0.1,c=hp_id_new)
# plt.colorbar()
# plt.show()

print('Making healpix')
# Create healpix map of sky
# nside is calculated as sqrt(npix/12) where npix is length of vectors
# order can be ring or nested, as determined by above plot and gorski healpix primer
hp_new = HEALPix(nside=64, order='nested',frame=Galactic())

print('Interpolating PSF')
# Get maximum angular extent of PSF
max_ang = np.int(np.floor(np.max(ang_dist)))

# Interpolate to highly-sampled binning
map_dist_lor = np.linspace(np.min(ang_dist),(max_ang),100000)
psf_interplor = interpolate.interp1d(ang_dist,brightness,kind='cubic')

print('Retrieving healpix pixels')
# Calculate omega and convert intensity to flux (W m^-2 sr^-1 -> W m^-2)
omega_new = ((hp_new.pixel_resolution.value/60)**2)*((np.pi/180)**2)
L_G_new_flux = L_G_new*omega_new # W m^-2

# Calculate Vega flux in W m^-2 at LORRI wavelength
lorri_nu = 3e8 / (607.62 * 1e-9) # LORRI pivot wavelength in Hz
vega_flux = 3050*(1e-26)*lorri_nu # Vega flux for LORRI band converted from Jy to W m^-2

# Make radius rings and add them up
def fun_4multiProc(argz):
    hp = argz[0] # Unpack the argz
    cent = argz[1]
    rad_bins_per = argz[2]
    rad_big = hp.cone_search_skycoord(cent,radius=rad_bins_per)
    print('Finished radius '+str(rad_bins_per))
    return rad_big

rad_bins = np.linspace(5,88,40,endpoint=True) # Make appropriate bins to plot
# Create psf_rad_bins
psf_rad_bins = psf_interplor(np.diff(rad_bins)/2+rad_bins[:-1]) # Get pts between bins for psf_rad_bins
psf_rad_bins_max = psf_rad_bins + (psf_rad_bins*0.1) # Maximum with 10% error
psf_rad_bins_min = psf_rad_bins - (psf_rad_bins*0.1) # Minimum with 10% error

#--- New Map ---
# Calculate flux in bins for new map
fun_vars = [None for i in range(0,rad_bins.size)]
for i in range(0,rad_bins.size):
    fun_vars[i] = (hp_new, cent, rad_bins[i]*u.degree) # Fill up function vars
# END FOR i
with Pool(7) as p:
    multi_done_new = p.map(fun_4multiProc, fun_vars) # Call the function in a multi-threaded scenario for speed
# END WITH
L_G_rad_new = np.zeros(rad_bins.size-1) # Preallocate
sumsum_new = 0
sumsum_pix_new = 0
for i in range(0,rad_bins.size-1,1):
    rad_lil = multi_done_new[i]
    rad_big = multi_done_new[i+1]
    if( np.sum(~np.isin(rad_big, rad_lil)) > 0 ):
        L_G_rad_new[i] = np.sum(L_G_new_flux[rad_big[~np.isin(rad_big, rad_lil)]])
    else:
        L_G_rad_new[i] = 0
    # END IF
    sumsum_new += L_G_rad_new[i]
    sumsum_pix_new += np.sum(~np.isin(rad_big, rad_lil))
# END FOR i
print('Total sum of new is: '+str(sumsum_new))
print('Total pix for new is: '+str(sumsum_pix_new))

# Center of rad bins
rad_bins_cent = np.diff(rad_bins)/2+rad_bins[:-1]

print('Saving data')
# Scale flux bins by psf beam in DN s^-1 pix^-1
L_G_scattered_new = L_G_rad_new/vega_flux*psf_rad_bins
L_G_scattered_new_fluxmaxerr = ((L_G_rad_new/vega_flux)+((L_G_rad_new/vega_flux)*0.0017))*psf_rad_bins
L_G_scattered_new_fluxminerr = ((L_G_rad_new/vega_flux)-((L_G_rad_new/vega_flux)*0.0017))*psf_rad_bins
L_G_scattered_new_psfmaxerr = L_G_rad_new/vega_flux*psf_rad_bins_max
L_G_scattered_new_psfminerr = L_G_rad_new/vega_flux*psf_rad_bins_min

# Convert from DN/s to nW m^-2 sr^-1
L_G_scattered_new_nW = L_G_scattered_new*sbconv
L_G_scattered_new_nW_psfmaxerr = L_G_scattered_new_psfmaxerr*sbconv
L_G_scattered_new_nW_psfminerr = L_G_scattered_new_psfminerr*sbconv
L_G_scattered_new_nW_fluxmaxerr = L_G_scattered_new_fluxmaxerr*sbconv
L_G_scattered_new_nW_fluxminerr = L_G_scattered_new_fluxminerr*sbconv

# Save important data to .mat file to read in with matlab
savemat(timestamp+'.mat', mdict={'scattered_nW': L_G_scattered_new_nW, 'rad_bins_cent' : rad_bins_cent,
                                 'psf_max_err' : L_G_scattered_new_nW_psfmaxerr, 'psf_min_err' : L_G_scattered_new_nW_psfminerr,
                                'flux_max_err' : L_G_scattered_new_nW_fluxmaxerr, 'flux_min_err' : L_G_scattered_new_nW_fluxminerr})

# # Plot summed flux bins
# plt.figure()
# # plt.scatter(rad_bins_cent,L_G_rad)
# plt.scatter(rad_bins_cent,L_G_rad_new)
# plt.xlabel(r'$\theta$ [$\degree$]')
# plt.ylabel('Bin-Summed Flux')
# plt.title(' Sum: '+'{0:.2f}'.format(np.sum(L_G_rad_new)))
# plt.show(block=False)

# # Plot centered rad bins with summed Vegas per bin
# plt.figure()
# # plt.scatter(rad_bins_cent,L_G_rad/vega_flux)
# plt.scatter(rad_bins_cent,L_G_rad_new/vega_flux)
# plt.xlabel(r'$\theta$ [$\degree$]')
# plt.ylabel('Bin-Summed Flux/Vega Flux')
# plt.title('Sum: '+'{0:.2f}'.format(np.sum(L_G_rad_new/vega_flux)))
# plt.show(block=False)

# Plot scattered intensity in nW and save in figs folder
# plt.figure()
# # plt.scatter(rad_bins_cent,L_G_scattered_nW)
# plt.xlabel(r'$\theta$ [$\degree$]')
# plt.ylabel(r'Bin-Summed Scattered Intensity [nW m$^{-2}$ sr$^{-1}$]')
# plt.scatter(rad_bins_cent,L_G_scattered_new_nW)
# plt.title('Sum: '+'{0:.2f}'.format(np.sum(L_G_scattered_new_nW)))
# plt.savefig('figs/'+str(timestamp)+'.png')
# plt.close()

print('done')