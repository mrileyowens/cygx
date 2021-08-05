# Authors: Loren Anderson, Riley Owens (GitHub: mrileyowens)

# This file creates a false-color image of Cygnus X in 8 and
# 24 micron bands overlaid by moment-zero RRL contours of the
# area and overplots HII regions from a catalog.

import matplotlib
from matplotlib import pyplot as plt

from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

from astropy.coordinates import SkyCoord
from regions import read_ds9,write_ds9,CirclePixelRegion,PixCoord,SkyRegion

import numpy as np

from reproject import reproject_interp

#Establishing directories
home='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU'
data=home+'/Data'
figs=home+'/Pictures and Plots'

def normalize(arr, vmin, vmax):
	nor = (arr - vmin) / (vmax - vmin)
	nor[np.where(nor<0.0)] = 0.0
	nor[np.where(nor>1.0)] = 1.0
	return nor

norm = (1.0, 200)

bl = "left"
ll = "bottom"

#coords = SkyCoord("79.5 0.8", frame="galactic", unit=(u.deg,u.deg))
#framesize = (u.Quantity(5.1,u.deg), u.Quantity(6.0,u.deg))

#Spitzer 8 and 24 micron maps and GBT moment-zero RRL map
fits_red=data+'/spitzer_M1_cygnus_2.4_24micron_reproj_cygnus_mosaic.fits'
fits_green=data+'/spitzer_I4_cygnus_2.4_8micron_reproj_cygnus_mosaic.fits'
fitsM0rrl=data+'/cygnus_m0_mosaic_-20_20.fits'

#Opening/creating data arrays and extracting WCS
red=fits.open(fits_red)
green=fits.open(fits_green)
blue=np.zeros(np.shape(red[0].data))
m0rrl=fits.open(fitsM0rrl)
#w_red = wcs.WCS(red[0].header)
#w_green=wcs.WCS(green[0].header)
w_m0rrl=wcs.WCS(m0rrl[0].header)

#cutout_red=Cutout2D(red[0].data, coords, framesize, wcs=w_red)
#cutout_green=Cutout2D(green[0].data, coords, framesize, wcs=w_green)
#cutout_blue=Cutout2D(np.zeros(np.shape(red[0].data)), coords, framesize, wcs=w_red)
#cutoutM0rrl=Cutout2D(m0rrl[0].data,coords,framesize,wcs=w_m0rrl)
#cutoutM0co=Cutout2D(m0co[0].data,coords,framesize,wcs=w_m0co)

#Normalizing false-color image arrays on a color scale
#wc = cutout_red.wcs
r = normalize(red[0].data, *norm)
g = normalize(green[0].data, *norm)
b = normalize(blue, *norm)

#Stacking normalized false-color image arrays
rgb=np.dstack([r,g,b])

#kernel=Gaussian2DKernel(x_stddev=7,y_stddev=7)
#cutoutM0coConv=convolve(cutoutM0co.data,kernel)

#Initializing figure
plt.close('all')
fig = plt.figure(figsize=(6,6))

#Adding false-color Spitzer maps to the figure
ax = fig.add_subplot(111, projection=w_m0rrl)
img = ax.imshow(rgb, origin="lower")

#Setting axes coordinates, labels, and tick format
lat=ax.coords['glat']
lat.set_axislabel('Galactic Latitude')
lat.set_major_formatter('d.d')
lon=ax.coords['glon']
lon.set_axislabel('Galactic Longitude')
lon.set_major_formatter('d.d')

#Overplotting moment-zero RRL contours
ax.contour(m0rrl[0].data,levels=[1,2,5,10],colors='white',alpha=0.5,linewidths=[0.5,1.0,1.5,2.0])
#ax.contour(m0rrl[0].data,levels=[1,2,4,8,16],colors='cyan',alpha=0.5,linewidths=[0.25,0.75,1.25,1.75,2.25],antialiased=True,nchunk=20)

#Overlaying HII regions in the area from existing catalog at http://astro.phys.wvu.edu/wise/
regFile=data+'/wise_matched_V2.3.reg'
regions=read_ds9(regFile)
for i in range(len(regions)):
    if (0<=i) and (i<=93):
        plotColor='red' #Known regions
    elif (94<=i) and (i<=124):
        plotColor='cyan' #Candidate regions
    elif (125<=i) and (i<=238):
        plotColor='yellow' #Radio-quiet regions
    else:
        plotColor='red' #More known regions
    pixel_region = regions[i].to_pixel(w_m0rrl)
    pixel_region.plot(color=plotColor)

fig.savefig(figs+'/cutRegions.png',dpi=200,bbox_inches='tight')
plt.show()
