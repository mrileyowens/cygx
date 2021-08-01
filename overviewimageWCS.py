import os
import sys
import shutil
import glob

import six

import matplotlib
#matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.ticker as ticker

from matplotlib import rc
#rc('text', usetex=True)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rcParams['font.size'] = 15
#plt.rcParams['text.latex.preamble'] = [r'\\usepackage[helvet]{sfmath}']
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'helvetica'

from astropy.io import fits
from astropy import wcs
#from astropy.table import Table, Column
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

from astropy.coordinates import SkyCoord
#from astropy.visualization import PercentileInterval
from astropy.visualization import AsinhStretch
from regions import read_ds9, write_ds9, RectangleSkyRegion,CirclePixelRegion,PixCoord,SkyRegion

import numpy as np

from reproject import reproject_interp

home='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU'
data=home+'/Data'
figs=home+'/Pictures and Plots'

def normalize(arr, vmin, vmax):
	nor = (arr - vmin) / (vmax - vmin)
	nor[np.where(nor<0.0)] = 0.0
	nor[np.where(nor>1.0)] = 1.0
	return nor

norm = (1.0, 200)
stretch = AsinhStretch(a=0.1)

bl = "left"
ll = "bottom"

sn = 2
sx = 2
sy = 1

coords = SkyCoord("79.5 0.8", frame="galactic", unit=(u.deg,u.deg)) # Galactic Centre
framesize = (u.Quantity(5.1,u.deg), u.Quantity(6.0,u.deg))

fits_red=data+'/spitzer_M1_cygnus_2.4_24micron_reproj_cygnus_mosaic.fits'
fits_green=data+'/spitzer_I4_cygnus_2.4_8micron_reproj_cygnus_mosaic.fits'
#fits_blue = "/home/loren/papers/gc_lobe/spitzer_gc_b.fits"
fitsM0rrl=data+'/cygnus_m0_mosaic_-20_20.fits'
#fitsM0co="C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/13CO_cube_CygX_m0_-20_20_reproj_spitzer_24micron.fits"

red=fits.open(fits_red)
green=fits.open(fits_green)
blue=np.zeros(np.shape(red[0].data))
#m0co=fits.open(fitsM0co)
m0rrl=fits.open(fitsM0rrl)
#blue = fits.open(fits_blue)
w_red = wcs.WCS(red[0].header)
w_green=wcs.WCS(green[0].header)
#w_m0co=wcs.WCS(m0co[0].header)
w_m0rrl=wcs.WCS(m0rrl[0].header)

#cutout_red=Cutout2D(red[0].data, coords, framesize, wcs=w_red)
#cutout_green=Cutout2D(green[0].data, coords, framesize, wcs=w_green)
#cutout_blue=Cutout2D(np.zeros(np.shape(red[0].data)), coords, framesize, wcs=w_red)
#cutoutM0rrl=Cutout2D(m0rrl[0].data,coords,framesize,wcs=w_m0rrl)
#cutoutM0co=Cutout2D(m0co[0].data,coords,framesize,wcs=w_m0co)

#wc = cutout_red.wcs
r = normalize(red[0].data, *norm)
g = normalize(green[0].data, *norm)
b = normalize(blue, *norm)

rgb=np.dstack([r,g,b])

#kernel=Gaussian2DKernel(x_stddev=7,y_stddev=7)
#cutoutM0coConv=convolve(cutoutM0co.data,kernel)



plt.close('all')
fig = plt.figure(figsize=(6,6))

ax = fig.add_subplot(111, projection=w_m0rrl)
img = ax.imshow(rgb, origin="lower")

#plt.show()

lat = ax.coords['glat']
lat.set_axislabel('Galactic Latitude')
lat.set_major_formatter('d.d')
lon = ax.coords['glon']
lon.set_axislabel('Galactic Longitude')
lon.set_major_formatter('d.d')

ax.contour(m0rrl[0].data,levels=[1,2,5,10],colors='white',alpha=0.5,linewidths=[0.5,1.0,1.5,2.0])
#ax.contour(m0rrl[0].data,levels=[1,2,4,8,16],colors='cyan',alpha=0.5,linewidths=[0.25,0.75,1.25,1.75,2.25],antialiased=True,nchunk=20)

#plt.show()

regFile = 'C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/wise_matched_V2.3.reg'
regions = read_ds9(regFile)
for i in range(len(regions)):
    if (0<=i) and (i<=93):
        plotColor='red'
    elif (94<=i) and (i<=124):
        plotColor='cyan'
    elif (125<=i) and (i<=238):
        plotColor='yellow'
    else:
        plotColor='red'
    print(i)
    print(plotColor)
    pixel_region = regions[i].to_pixel(w_m0rrl)
    pixel_region.plot(color=plotColor)  # plot options passed to matplotlib



fig.savefig(figs+'/cutRegions.png',dpi=200,bbox_inches='tight')
plt.show()
#fig.savefig("/home/loren/papers/gc_lobe/spitzer_rgb.pdf",dpi=1000,pad_inches=0.2,bbox_inches='tight')
