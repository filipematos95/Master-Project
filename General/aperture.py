from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import CircularAperture
import numpy as np
from photutils import aperture_photometry
from astropy.io import fits
from photutils import CircularAnnulus
from photutils import Background2D, MedianBackground
from astropy.stats import SigmaClip
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize,LinearStretch,ZScaleInterval)

"""

This is a code performs aperture photometry on a single star.
This code is not a standalone but a module. Should be used with aperture_photometry.py

"""

def aper_phot(img_data,positions, r = 10., r_in = 14., r_out = 18,bkg_sub=False,plot=False):
"""
:params: r: Aperture Radius
:params: r_in: annulus aperture inside radius
:params: r_out: annulus aperture outside radius
:params: bkg_sub: True if background subtraction is needed
:params: plot: True to plot

:out: phot_table: Table with the values of the aperture photometry
"""

	#Background subtraction
	if bkg_sub == True:
		sigma_clip = SigmaClip(sigma=3., iters=10)
		bkg_estimator = MedianBackground()
		bkg = Background2D(img_data, (50, 50), filter_size=(3, 3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
		data_sub = img_data - bkg.background
	else:
		data_sub = img_data

	#Aperture Photometry using a circular aperture and a circular annulus
	apertures = CircularAperture(positions, r=r)
	annulus_apertures = CircularAnnulus(positions,r_in = r_in,r_out = r_out)
	apers = [apertures,annulus_apertures]
	phot_table = aperture_photometry(data_sub, apers)
	bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
	bkg_sum = bkg_mean * apertures.area()
	final_sum = phot_table['aperture_sum_0'] - bkg_sum
	phot_table['residual_aperture_sum'] = final_sum
	positions = np.array(positions)

	if plot == True:
		#Ploting 
		norm = ImageNormalize(data_sub, interval=ZScaleInterval(),stretch=LinearStretch())
		plt.imshow(data_sub, cmap='Greys', origin='lower',norm=norm)
		apertures.plot(color='blue', lw=1.5, alpha=0.5)
		annulus_apertures.plot(color='green',lw=1.5,alpha=0.5)
		plt.plot(positions[:,0],positions[:,1], ls='none', color = 'red', marker='.', ms=10, lw=1.5)
		plt.xlim(0, data_sub.shape[1]-1)
		plt.ylim(0, data_sub.shape[0]-1)
		plt.show()	

	return phot_table