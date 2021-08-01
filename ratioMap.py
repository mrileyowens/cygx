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

regFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/rrlRegion.reg'
file1='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/cygnus_m0_mosaic_-20_20_reproj_CGPS_1420_MHz.fits'
file2='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/CGPS_CygXc_1420_MHz_cnv2p0.fits'

hdul1=fits.open(file1)
hdul2=fits.open(file2)

wcs=wcs.WCS(hdul1[0].header)

contData=hdul2[0].data+6.0
contData=contData*0.0513

ratio=np.divide(hdul1[0].data,np.squeeze(contData))
ratio[np.where(hdul1[0].data<1.0)]=np.nan
#ratio=ratio[:,~np.isnan(hdul1[0].data).all(axis=0)]
#ratio=ratio[~np.isnan(hdul1[0].data).all(axis=1)]

#contourData=hdul1[0].data[:,~np.isnan(hdul1[0].data).all(axis=0)]
#contourData=contourData[~np.isnan(contourData).all(axis=1)]

#ratioNorm=normalize(ratio,*norm)
#ratioNorm[:,~np.isnan(ratioNorm).all(axis=1)]

plt.close('all')
fig=plt.figure()

ax=fig.add_subplot(111,projection=wcs)
img=ax.imshow(ratio,origin='lower',cmap='inferno',vmin=0.0,vmax=5.0)

#ax.contour(hdul1[0].data,levels=[1.0],colors='black',alpha=0.5,linewidths=[1.0])

fig.colorbar(img,ax=ax,label='Ratio (RRL / Continuum)')

lat=ax.coords['glat']
lat.set_axislabel('Galactic Latitude (deg.)')
lat.set_major_formatter('d.d')
lon=ax.coords['glon']
lon.set_axislabel('Galactic Longitude (deg.)')
lon.set_major_formatter('d.d')

region=read_ds9(regFile)
pixelRegion=region[1].to_pixel(wcs)
pixelRegion.plot(color='black',linestyle='dashed')

plt.tight_layout()

plt.show()

plt.close('all')

fig2=plt.figure()

ax2=fig2.add_subplot(111,projection=wcs)

eTemp=(7103.3*(5.7578**1.1)*(np.divide(np.squeeze(contData),hdul1[0].data))*((1.07)**(-1.0)))**0.87
eTemp[np.where(hdul1[0].data<1.0)]=np.nan

img2=ax2.imshow(eTemp,origin='lower',cmap='inferno',vmin=0.0,vmax=10000.0)
fig2.colorbar(img2,ax=ax2,label='Electron Temperature (K)')

lat2=ax2.coords['glat']
lat2.set_axislabel('Galactic Latitude (deg.)')
lat2.set_major_formatter('d.d')
lon2=ax2.coords['glon']
lon2.set_axislabel('Galactic Longtiude (deg.)')
lon2.set_major_formatter('d.d')

region=read_ds9(regFile)
pixelRegion=region[1].to_pixel(wcs)
pixelRegion.plot(color='black',linestyle='dashed')

plt.tight_layout()

plt.show()

eTemp=eTemp[~np.isnan(eTemp)]
counts,bins=np.histogram(eTemp,bins=200)
plt.hist(bins[:-1],bins,weights=counts)
plt.xlabel('Electron Temperature (K)')
plt.ylabel('Pixel Count')
plt.title('Electron Temperature Histogram Distribution')
plt.xlim(xmin=0,xmax=20000)

plt.show()
