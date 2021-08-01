from spectral_cube import SpectralCube
import astropy.units as u

#Creating spectral cube from 13CO emission datacube
datacubeFile='C:/Users/15136/OneDrive - University of Cincinnati/Documents/Research/WVU REU/Data/13CO_cube_CygX.fits'
cube=SpectralCube.read(datacubeFile)

#Assigning spectral unit to cube
cube=cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')

#Selecting a subset of the cube
cubeSub=cube.spectral_slab(-20.0*u.km/u.s,20.0*u.km/u.s)

#Creating moment-one map from subcube
cubeSubM1=cubeSub.moment(order=1)

#Saving moment-one map
cubeSubM1.write('13CO_cube_CygX_m1_-20_20.fits',overwrite=True)
