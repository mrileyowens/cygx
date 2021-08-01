# Author: Riley Owens (GitHub: mrileyowens)

# This file creates a colorscaled moment-zero
# RRL map with a colorbar overlaid by square
# tiles showing where the C-band RRL cubes
# have been processed by GaussPy+.

from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from reproject import reproject_interp

import matplotlib
from matplotlib import pyplot as plt

import numpy as np

from regions import read_ds9, write_ds9, RectangleSkyRegion,CirclePixelRegion,PixCoord,SkyRegion

# Establishing directories and filepaths
home='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU'
data=home+'/Data'
figs=home+'/Pictures and Plots'

fitsRRLm0=data+'/cygnus_m0_mosaic_-20_20.fits'
fitsRRLm1=data+'/cygnus_m1_mosaic_-20_20.fits'

# Setting HDU lists
hdulM0=fits.open(fitsRRLm0)
hdulM1=fits.open(fitsRRLm1)

# Extracting moment-zero/one RRL data
m0data=hdulM0[0].data
m1data=hdulM1[0].data

# Reassigning low-signal areas in the moment-zero
# data as NaNs in the moment-one data
#m1data[np.where(m0data<0.5)]=np.nan

# Extracting the WCS to be used
wcs=wcs.WCS(hdulM1[0].header)

# Initializing figure
plt.close('all')
fig=plt.figure()

# Adding subplot of moment-zero RRL data on colormap between 0-5 K km/s
ax=fig.add_subplot(111,projection=wcs)
img=ax.imshow(m0data,origin='lower',cmap='gray',vmin=0.0,vmax=5.0,aspect='equal')

# Setting coordinates, labels, and tick format of axes
lat=ax.coords['glat']
lat.set_axislabel('Galactic Latitude (deg.)')
lat.set_major_formatter('d.d')
lon=ax.coords['glon']
lon.set_axislabel('Galactic Longitude (deg.)')
lon.set_major_formatter('d.d')

# Adding colorbar with same height as subplots
cbar=plt.colorbar(img,ax=ax,fraction=0.047*(m0data.shape[0]/m0data.shape[1]))
cbar.set_label('Intensity (K km s$^{-1}$)',rotation=270,labelpad=18.0)

# Overplotting rectangular regions where RRL data has been decomposed by GaussPy+
regFile = 'C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/agdRegions.reg'
regions = read_ds9(regFile)
for i in range(len(regions)):
    pixel_region = regions[i].to_pixel(wcs)
    pixel_region.plot(color='yellow')

#fig.colorbar(img,ax=ax,fraction=0.047*(m1data.shape[0]/m1data.shape[1]))

#plt.tight_layout()

plt.savefig(figs+'/m0rrl.png',dpi=1000,bbox_inches='tight',overwrite=True)

plt.show()
