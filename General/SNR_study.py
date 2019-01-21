import numpy as np
from photutils import IRAFStarFinder
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils import datasets
import matplotlib.pyplot as plt
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize,LinearStretch,ZScaleInterval)
from photutils import CircularAperture
from astropy.io import fits
from aperture import aper_phot
from photutils import find_peaks
from photutils.datasets import make_100gaussians_image
from photutils import Background2D, MedianBackground
from astropy.stats import SigmaClip
from photutils.psf import IterativelySubtractedPSFPhotometry
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from mpl_toolkits.mplot3d import Axes3D
import time
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

DIR = #Directory with the files to load
DIR_TO_SAVE = #Directory to save the files

name = 'NAME OF THE TARGET'

with open(DIR+'list.txt', 'r') as fp: #List with the filenames to analyse
	files = fp.read().splitlines()

files_not_use = []
bad_files = np.load('bad_files_'+name+'.npy') #File with the files to skip (due to bad quality)

def star_fits (filename,sort):
	hdulist = fits.open(filename,ignore_missing_end=True)
	data = hdulist[1].data
	data = data[data[sort].argsort()]

	return data

datas = star_fits('filename.fit','Jmag')
positions_stars = []
stars = []

for dat in datas[0:5]:
	ra = dat['RAJ2000'] *u.deg
	dec = dat['DEJ2000'] * u.deg
	positions_stars.append(SkyCoord(ra,dec,frame = 'icrs'))
	stars.append(dat['_2MASS'])

plot = True #True to actiave interactive plotting
appertures = []
appertures_error = []
date = []
files = np.array(files)
arg_del = np.in1d(files,bad_files)
arg_del = np.invert(arg_del)
files = files[arg_del]
SNR_r = []
noise_r = []
appertures_r = []
rn = np.random.randint(300,size=10)
just_to_know = 0
time0 = time.time()
SNR_best = []

stars = stars[0]

SNR_r = []
noise_r = []
appertures_r = []
just_to_know = 0
SNR_best = []
appertures = []
appertures_error = []
date = []
position = positions_stars[0]
best_r = {}

for filename in files: 
	print just_to_know
	just_to_know += 1
	try:
		hdulist = fits.open(DIR+filename,ignore_missing_end=True)
		data = hdulist['SCI'].data
		header = hdulist['SCI'].header
		rdnoise = header['RDNOISE']
		hdulist.close()
	except:
		print filename
		continue

	print "Working: " + filename

 	size = u.Quantity((50.0, 50.0), u.arcsec)
	cutout = Cutout2D(data, position, size, wcs=WCS(header))
	image_cutted = cutout.data
	
	mjd = header['MJD-OBS']

	mean, median, std = sigma_clipped_stats(image_cutted,sigma = 5.0, iters = 10) 
	threshold = median + (2.0 * std) #Here the user should insert the threshold to find the centroid of the star. Default is 2 sigma from the median
	fwhm = 10.0 # The user should insert an approximate FWHM.
	sigma_psf = float(fwhm)/gaussian_sigma_to_fwhm

	#Finding the centroid of the star
	daofind = DAOStarFinder(fwhm=sigma_psf*gaussian_sigma_to_fwhm, threshold=2*threshold,sky = median) 

	sources = daofind(image_cutted-median)
	argsources = sources.argsort('flux') 
	argsources = np.flip(argsources,0)

	positions = []
	SNR = []
	rs = []

	appertures_temp = []
	appertures_error_temp = []
	date_temp = []

	#Computes the SNR for apertures sizes from 3 to 20 pixels. 
	for r in range(3,20):
		if True:
			if len(sources['xcentroid']) > 0:
				if name == 'WASP-80': #Filter to choose the WASP-80 instead of its binary companion
					if sources['xcentroid'][0] < 75:
						source_x = sources['xcentroid'][0]
						source_y = sources['ycentroid'][0]
						positions.append([sources['xcentroid'][0], sources['ycentroid'][0]])
					else:
						source_x = sources['xcentroid'][1]
						source_y = sources['ycentroid'][1]
						positions.append([sources['xcentroid'][1], sources['ycentroid'][1]])
				else:
					source_x = sources['xcentroid'][0]
					source_y = sources['ycentroid'][0]
					positions.append([sources['xcentroid'][0], sources['ycentroid'][0]])

				aper =  aper_phot(image_cutted,positions,r = r, r_in = 40, r_out = 45 ,bkg_sub = False,plot=False)
				if np.isnan( np.sqrt(-aper['residual_aperture_sum'][0] + aper['aperture_sum_0'][0] )) == False:
					N_star = aper['residual_aperture_sum'][0]
					N_bkg = (-aper['residual_aperture_sum'][0] + aper['aperture_sum_0'][0]) 
					N_rdout =  rdnoise*rdnoise*np.pi*r*r
					noise2 = np.sqrt(N_star + N_bkg + N_rdout)
					SNR.append(aper['residual_aperture_sum'][0] / noise2)
					appertures_temp.append(aper['residual_aperture_sum'][0])
					appertures_error_temp.append(np.sqrt(aper['residual_aperture_sum'][0] + rdnoise**2))
					date_temp.append(mjd)
					rs.append(r)
			else:
				print 'No stars found in the file: ' + filename
			
		else:
			print 'No stars found in the file(except): ' + filename
	
	if len(SNR) > 0:	
		best = 	np.argmax(SNR)
		SNR_best.append(rs[best])
		
		appertures.append(appertures_temp[best])
		appertures_error.append(appertures_error_temp[best])
		date.append(date_temp[best])

		SNR = np.array(SNR)
		SNR_mean = np.mean(SNR)
		SNR_median = np.mean(SNR)
		SNR_r.append([r,SNR_mean,SNR_median])
		
		best_r[filename] = rs[best]

		print rs[best]

np.save(DIR_TO_SAVE+'SNR_radius_dic_'+name+'.npy',best_r)

print "Time elapsed: " + str((time.time() - time0))

#Plotting and savng the results
appertures = np.array(appertures)
appertures_error = np.array(appertures_error)
date = np.array(date)
plt.errorbar(date,appertures,yerr = appertures_error,fmt='o')
plt.title('Flux ' + name)
plt.xlabel('Date [MJD]')
plt.ylabel('Flux [Counts]')
plt.savefig(DIR_TO_SAVE+'Flux_'+name+'_'+str(r)+'_out_error.pdf')
plt.close()
np.save(DIR_TO_SAVE+'Flux_'+name+'_'+str(r)+'_out_error.npy',(date,appertures,appertures_error))

print "Median: " + str(np.mean(SNR_best))
print "Mean: " + str(np.median(SNR_best))

bins = np.linspace(2.5,15.5,14)
np.save(DIR_TO_SAVE+'SNR_radius_'+name+'.npy',SNR_best)
plt.hist(SNR_best,bins = bins,width = 1)
plt.title('SNR distribution')
plt.xlabel('Radius [pix]')
plt.ylabel('SNR')
plt.savefig(DIR_TO_SAVE+'Histogram_SNR.pdf')
plt.show()
plt.clf()

'''
SNR_r = np.array(SNR_r)
np.save(DIR_TO_SAVE+'SNR_r_'+name+'_noise_v2.npy',SNR_r)
plt.title("SNR: MEAN")
plt.plot(SNR_r[:,0],SNR_r[:,1],color='blue')
plt.xlabel("Radius [pix]")
plt.ylabel("SNR")
#plt.savefig(DIR_TO_SAVE+'SNR_mean_'+name+'_noise_v2.pdf')
#plt.show()
plt.clf()

plt.plot(SNR_r[:,0],SNR_r[:,2],color = 'red')
plt.xlabel("Radius [pix]")
plt.ylabel("SNR")
plt.title("SNR: MEDIAN")
#plt.savefig(DIR_TO_SAVE+'SNR_median_'+name+'_noise_v2.pdf')
plt.clf()
'''

