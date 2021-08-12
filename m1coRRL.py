# Author: Riley Owens (GitHub: mrileyowens)

# This file creates a side-by-side figure
# of moment-one 13CO and RRL data in
# Cygnus X with corresponding areas of
# low signal in the moment-zero data masked.

from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from reproject import reproject_interp

import matplotlib
from matplotlib import pyplot as plt

import numpy as np

# Establishing directories
home='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU'
data=home+'/data'

fitsRRLm0=data+'/cygnus_m0_mosaic_-20_20.fits'
fitsRRLm1=data+'/cygnus_m1_mosaic_-20_20.fits'
fitsCOm0=data+'/13CO_cube_CygX_m0_-20_20.fits'
fitsCOm1=data+'/13CO_cube_CygX_m1_-20_20.fits'

# Setting HDU lists
hdul1=fits.open(fitsRRLm0)
hdul2=fits.open(fitsRRLm1)
hdul3=fits.open(fitsCOm0)
hdul4=fits.open(fitsCOm1)

# Extracting moment-zero/one 13CO data
coM0data=hdul3[0].data
coM1data=hdul4[0].data

# Creating a dummy array to be used as a mask
coM1mask=np.ones(np.shape(coM0data))
rrlM1mask=np.ones(np.shape(coM0data))

# Reprojecting/interpolating RRL data to the WCS of the 13CO data
rrlM0dataReproj,footprintM0=reproject_interp(hdul1[0],hdul4[0].header)
rrlM1dataReproj,footprintM1=reproject_interp(hdul2[0],hdul4[0].header)

# Reassigning low-signal areas in the moment-zero RRL
# map as NaNs in the moment-one RRL map
rrlM1dataReproj[np.where(rrlM0dataReproj<0.5)]=np.nan

# Reassigning low-signal areas in the moment-zero 13CO
# map as NaNs in the moment-one 13CO map
coM1data[np.where(coM0data<0.5)]=np.nan

# Repeating the low-signal masking for the masks
coM1mask[np.where(coM0data<0.5)]=0.0
rrlM1mask[np.where(rrlM0dataReproj<0.5)]=0.0

# Creating HDU from masked moment-one data
hdu2=fits.PrimaryHDU(rrlM1dataReproj,header=hdul4[0].header)
hdu4=fits.PrimaryHDU(coM1data,header=hdul4[0].header)

# Creating HDU from low-signal masks
hduCOmask=fits.PrimaryHDU(coM1mask,header=hdul4[0].header)
hduRRLmask=fits.PrimaryHDU(rrlM1mask,header=hdul4[0].header)

# Creating HDU lists from masked moment-one HDUs
hdul2reduced=fits.HDUList([hdu2])
hdul4reduced=fits.HDUList([hdu4])

# Creating HDU lists from mask HDUs
hdulCOmask=fits.HDUList([hduCOmask])
hdulRRLmask=fits.HDUList([hduRRLmask])

# Saving the masked moment-one maps
#hdul2reduced.writeto('cygnus_m1_mosaic_20_20_reproj_13CO_cube_CygX_reduced.fits',overwrite=True)
#hdul4reduced.writeto('13CO_cube_CygX_m1_-20_20_reduced.fits',overwrite=True)

# Saving the masks
#hdulRRLmask.writeto('cygnus_m1_mosaic_-20_20_reproj_13CO_cube_CygX_reduced_mask.fits',overwrite=True)
#hdulCOmask.writeto('13CO_cube_CygX_m1_-20_20_reduced_mask.fits')

# Extracting the WCS we will use in the figure
wcs1=wcs.WCS(hdul4[0].header)

# Initializing figure
plt.close('all')
fig=plt.figure()

# Adding two subplots to the figure placed next to each other
# with the same WCS and axes
ax1=fig.add_subplot(121,projection=wcs1)
ax2=fig.add_subplot(122,projection=wcs1,sharex=ax1,sharey=ax1)

# Plotting the moment-one RRL and 13CO data on a colormap between +/- 10 km/s
img1=ax1.imshow(coM1data,origin='lower',cmap='inferno',vmin=-10.0,vmax=10.0,aspect='equal')
img2=ax2.imshow(rrlM1dataReproj,origin='lower',cmap='inferno',vmin=-10.0,vmax=10.0,aspect='equal')

# Setting the axes coordinates, labels, and
# tick formats for both subplots
lat1=ax1.coords['ra']
lat1.set_axislabel('Right Ascension (deg.)')
lat1.set_major_formatter('d.d')
lon1=ax1.coords['dec']
lon1.set_axislabel('Declination (deg.)')
lon1.set_major_formatter('d.d')

lat2=ax2.coords['ra']
lat2.set_axislabel('Right Ascension (deg.)')
lat2.set_major_formatter('d.d')
lon2=ax2.coords['dec']
lon2.set_axislabel('')
lon2.set_major_formatter('d.d')
lon2.set_ticklabel_visible(False)

#overlay=ax2.get_coords_overlay('galactic')
#overlay.grid(color='lime',ls='dashed')
#overlay[0].set_axislabel('Galactic Longitude')
#overlay[1].set_axislabel('Galactic Latitude')

# Creating a colorbar the same height as the subplots
fig.colorbar(img2,ax=ax2,fraction=0.047*(rrlM1dataReproj.shape[0]/rrlM1dataReproj.shape[1]))

plt.tight_layout()

plt.show()
