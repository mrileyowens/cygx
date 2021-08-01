# Author: Riley Owens (GitHub: mrileyowens)

# This file takes one map and reprojects and
# interpolates it to a different map's WCS.

from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import numpy as np

#Establishing directories
home='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU'
data=home+'/Data'
figs=home+'/Pictures and Plots'

# The image to be reprojected and interpolated
fitsFile1=data+'/spitzer_M1_cygnus_2.4_24micron.fits'

# The image to reproject to and interpolate across
fitsFile2=data+'/cygnus_m0_mosaic_-20_20.fits'

#Setting HDU lists
hdul1=fits.open(fitsFile1)
hdul2=fits.open(fitsFile2)

#hdul2WCS=WCS(hdul2[0].header)

#hdul2WCS=WCS.dropaxis(hdul2WCS,-1)
#hdul2WCS=WCS.dropaxis(hdul2WCS,-1)

#hdul2WCSheader=hdul2WCS.to_header()

#Performing reprojection and interpolation
reprojArray,footprint=reproject_interp(hdul1[0],hdul2[0].header)

#Saving new reprojected/interpolated file
hdu=fits.PrimaryHDU(reprojArray,header=hdul2[0].header)
hdul3=fits.HDUList([hdu])
hdul3.writeto(data+'/spitzer_M1_cygnus_2.4_24micron_reproj_cygnus_mosaic.fits',overwrite=True)
