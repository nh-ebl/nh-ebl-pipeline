import sys
import os
from astropy.io import fits
sys.path.append('../')
from gdpyc import DustMap
import numpy as np
from astropy.wcs import WCS

#this is a reference image
# hdul = fits.open('regist_cruise_20081016_008642_lor_0086422037_0x633_sci_1.fit')
#<read in the fits image name from a text file made by matlab>
currentDir = os.path.dirname(os.path.realpath(__file__)) #get current working directory (calling from cmd line python has directory issues)
os.chdir(currentDir) #change dir to the needed directory
fileID = open("imagefile.txt","r") #open the file
imagefile = fileID.readlines()[0] #get the imagefile path
fileID.close() #close the file
hdul = fits.open(imagefile) #use the imagefile path

# Get header info
ref_head = hdul[0].header
timestamp = ref_head['SPCUTCJD'][3:]

# Setting up WCS
w = WCS(hdul[0].header)

# Make grid of coords
xx, yy = np.meshgrid(np.arange(0,ref_head['NAXIS1'],1),np.arange(0,ref_head['NAXIS2'],1))
sky = w.pixel_to_world(xx, yy)

# Calculate extinction map
list_filters = DustMap._parse_filters('Landolt_R') # Accepts a string (one filter) or list of strings (multiple filters)
ebv = DustMap.ebv(sky, dustmap='SFD', hires=True) # Calculate reddening
aebv = DustMap._ebv_to_ext(list_filters) # Get coefficient for filter (using Rv = 3.1 - can't change)
ext_map = np.squeeze(np.matmul(np.expand_dims(ebv, axis=-1), np.expand_dims(aebv, axis=0))) # Multiplies ebv*aebv but can support multiple filters (not just 1)
#ext_map = DustMap.extinction(sky, dustmap='SFD', filters='Landolt_R', hires=True) # Doesn't work due to matrix multiplication error

# import matplotlib.pyplot as plt
# plt.figure()
# im = plt.imshow(ext_map)
# plt.colorbar(im,label=r'Extinction')
# plt.xlabel(r'X Pix')
# plt.ylabel(r'Y Pix]')
# plt.show()

# Save to fits file
hdu1 = fits.PrimaryHDU(ext_map)
hdul1 = fits.HDUList([hdu1])
hdul1.writeto('ext_' + timestamp + '.fits',overwrite=True)

print('done')