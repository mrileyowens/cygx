# Author: Riley Owens (GitHub: mrileyowens)

# This file creates a moment-one map of 13CO emission data
# from Cygnus X overlaid by moment-one RRL contours in the
# same region, following a red/blueshift convention.

import matplotlib
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import wcs

import numpy as np

# Establishing directories
fitsCO='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/13CO_cube_CygX_m1_-20_20_reduced.fits'
fitsRRL='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/cygnus_m1_mosaic_-20_20_reproj_13CO_cube_CygX_reduced.fits'

# Setting the HDU lists of the maps
hdulCO=fits.open(fitsCO)
hdulRRL=fits.open(fitsRRL)

# Extracting the WCS from the 13CO map
wcs=wcs.WCS(hdulCO[0].header)

# Extracting moment-one RRL data and droppning rows/columns where
# the moment-one 13CO map is not defined
rrlData=hdulRRL[0].data[:,~np.all(np.isnan(hdulCO[0].data),axis=0)]
rrlData=rrlData[~np.all(np.isnan(hdulCO[0].data),axis=1)]

# Extracting moment-one 13CO data and masking rows/columns where
# it is not defined
coData=hdulCO[0].data[:,~np.all(np.isnan(hdulCO[0].data),axis=0)]
coData=coData[~np.all(np.isnan(coData),axis=1)]

# Initializing the figure
plt.close('all')
fig=plt.figure()
ax=fig.add_subplot(111,projection=wcs)

# Plotting moment-one 13CO map with a colormap between +/-10 km/s
img=ax.imshow(coData,origin='lower',cmap='spring',vmin=-10.0,vmax=10.0)

# Setting coordinates, labels, and tick
# format of axes
lat=ax.coords['ra']
lat.set_axislabel('Right Ascension (deg.)')
lat.set_major_formatter('d.d')
lon=ax.coords['dec']
lon.set_axislabel('Declination (deg.)')
lon.set_major_formatter('d.d')

#Adding contours of the moment-one RRL data using red/blueshift convention
ax.contour(rrlData,levels=[5.0,10.0],colors='cyan',alpha=1.0,linewidths=[0.5,1.0])
ax.contour(rrlData,levels=[-10.0,-5.0],colors='red',alpha=1.0,linewidths=[0.5,1.0])
fig.colorbar(img,ax=ax)

plt.show()
