# Author: Riley Owens (GitHub: mrileyowens)

# This file creates a ratio map between
# moment-zero RRL and radio continuum Data
# of Cygnus X, which is used to create a
# map of the electron temperature in the
# region.

from astropy.io import fits
from astropy import units as u
from astropy import wcs

from regions import read_ds9, write_ds9, RectangleSkyRegion,CirclePixelRegion,PixCoord,SkyRegion

import numpy as np

from matplotlib import pyplot as plt

def normalize(arr,vmin,vmax):
    nor=(arr-vmin)/(vmax-vmin)
    nor[np.where(nor<0.0)]=0.0
    nor[np.where(nor>1.0)]=1.0
    return nor

norm=(1.0,5.0)

# Establishing filepaths
regFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/rrlRegion.reg'
file1='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/cygnus_m0_mosaic_-20_20_reproj_CGPS_1420_MHz.fits'
file2='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/CGPS_CygXc_1420_MHz_cnv2p0.fits'

# Setting HDU lists
hdul1=fits.open(file1)
hdul2=fits.open(file2)

# Extracting WCS
wcs=wcs.WCS(hdul1[0].header)

# Adjusting background continuum and scaling
# the continuum data
contData=hdul2[0].data+6.0
contData=contData*0.0513

# Creating ratio map between continuum data and RRL data
# and reassigning low-signal areas in the moment-zero RRL
# data as NaNs in the ratio map
ratio=np.divide(hdul1[0].data,np.squeeze(contData))
ratio[np.where(hdul1[0].data<1.0)]=np.nan
#ratio=ratio[:,~np.isnan(hdul1[0].data).all(axis=0)]
#ratio=ratio[~np.isnan(hdul1[0].data).all(axis=1)]

#contourData=hdul1[0].data[:,~np.isnan(hdul1[0].data).all(axis=0)]
#contourData=contourData[~np.isnan(contourData).all(axis=1)]

#ratioNorm=normalize(ratio,*norm)
#ratioNorm[:,~np.isnan(ratioNorm).all(axis=1)]

# Initializing figure
plt.close('all')
fig=plt.figure()

# Plotting colorscaled ratio map
ax=fig.add_subplot(111,projection=wcs)
img=ax.imshow(ratio,origin='lower',cmap='inferno',vmin=0.0,vmax=5.0)

#ax.contour(hdul1[0].data,levels=[1.0],colors='black',alpha=0.5,linewidths=[1.0])

# Adding colorbar to the figure
fig.colorbar(img,ax=ax,label='Ratio (RRL / Continuum)')

# Setting coordinates, labels, and tick format to axes
lat=ax.coords['glat']
lat.set_axislabel('Galactic Latitude (deg.)')
lat.set_major_formatter('d.d')
lon=ax.coords['glon']
lon.set_axislabel('Galactic Longitude (deg.)')
lon.set_major_formatter('d.d')

# Overplotting the area of the moment-zero RRL data
region=read_ds9(regFile)
pixelRegion=region[1].to_pixel(wcs)
pixelRegion.plot(color='black',linestyle='dashed')

plt.tight_layout()
plt.show()

# Initializing figure
plt.close('all')
fig2=plt.figure()

# Adding subplot
ax2=fig2.add_subplot(111,projection=wcs)

# Creating an electron temperature map from the ratio map based on Eq. 1 from Quireza et al. 2006
# (https://iopscience.iop.org/article/10.1086/508803/pdf)
eTemp=(7103.3*(5.7578**1.1)*(np.divide(np.squeeze(contData),hdul1[0].data))*((1.07)**(-1.0)))**0.87
eTemp[np.where(hdul1[0].data<1.0)]=np.nan

# Plotting the colorscaled electron temperature map
img2=ax2.imshow(eTemp,origin='lower',cmap='inferno',vmin=0.0,vmax=10000.0)
fig2.colorbar(img2,ax=ax2,label='Electron Temperature (K)')

# Setting coordinates, labels, and tick format of axes
lat2=ax2.coords['glat']
lat2.set_axislabel('Galactic Latitude (deg.)')
lat2.set_major_formatter('d.d')
lon2=ax2.coords['glon']
lon2.set_axislabel('Galactic Longtiude (deg.)')
lon2.set_major_formatter('d.d')

# Overplotting the area of the moment-zero RRL data
region=read_ds9(regFile)
pixelRegion=region[1].to_pixel(wcs)
pixelRegion.plot(color='black',linestyle='dashed')

plt.tight_layout()
plt.show()

# Creating a histogram showing distribution of electron
# temperature in the map
eTemp=eTemp[~np.isnan(eTemp)]
counts,bins=np.histogram(eTemp,bins=200)
plt.hist(bins[:-1],bins,weights=counts)
plt.xlabel('Electron Temperature (K)')
plt.ylabel('Pixel Count')
plt.title('Electron Temperature Histogram Distribution')
plt.xlim(xmin=0,xmax=20000)

plt.show()
