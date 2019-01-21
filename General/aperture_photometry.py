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
#from mpl_toolkits.mplot3d import Axes3D
import time
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import os
import matplotlib.colors as colors

DIR = #Directory with the files to load (frames)
DIR_TO_SAVE_FLUX = #Directory to save te results from aperture photometry for all stars
DIR_TO_SAVE_BACKGROUND = #Directory to save the background test results for the entire frame 
DIR_TO_SNR = #Directory with the list of optimal aperture sizes for the target
stars_file = 'filename.fit' #Fits file with the list of stars 
files_list = 'list.txt' #Text file with the names of the frames files

bkg_study = False #True to preform background tests
fluxes_study = True #True to preform aperture photometry
bad_files = False #True to produce the bad_files. This is a manual task where the user views all the imeages and clicks on the bad ones generating a file with the names of the images to skip.
plot = False #Interactive plotting of the frames without apertures. Mandatory to produce the bad files list
plot_aper = False #Interactive plotting of the appertures. Usefull to check the performance of the code.
bkg_sub = False #Performs background subtraction instead of using an anullus aperture. Usefull for crowded fields. 
r = 14 #Default aperture radius
days = np.array([])
global target_name 
target_name = #Name of the target

#Check directories
if not os.path.exists(DIR_TO_SAVE_FLUX):
    os.makedirs(DIR_TO_SAVE_FLUX)

if not os.path.exists(DIR_TO_SAVE_BACKGROUND):
    os.makedirs(DIR_TO_SAVE_BACKGROUND)

files_not_use = []

#Button click event function
def onclick(event):
    print('Bad Picture')
    files_not_use.append(filename)
    print filename
    return 1 

#Opens the fits file with the star information
def star_fits (filename,sort):
	hdulist = fits.open(filename,ignore_missing_end=True)
	data = hdulist[1].data
	data = data[data[sort].argsort()]

	return data

#Selects the images that are not in the bad files list
def open_files(files_list,name = target_name):
	with open(DIR+files_list, 'r') as fp:
		files = fp.read().splitlines()

	if not os.path.exists('bad_files_'+name+'.npy'):
		np.save('bad_files_'+name+'.npy',np.array([]))
	bad_files = np.load('bad_files_'+name+'.npy')

	files = np.array(files)
	arg_del = np.in1d(files,bad_files)
	arg_del = np.invert(arg_del)
	files = files[arg_del]

	return files

def choose_stars(starsfile,n_stars = 5):
	'''
	Creates an array with the position of all the stars and their names. Selects which stars to use by taking the n_stars with a lower magnitude after sorting. 

	:params: starsfile: name of the file with all the stars
	:params: n_stars: number of stars to analyse. 
	'''
	datas = star_fits(starsfile,'Bmag')
	datas = datas[np.where(datas['Bmag'] > 0.0)]
	positions_stars = []
	stars = []
	for dat in datas[:n_stars]:
		ra = dat['RAJ2000'] *u.deg
		dec = dat['DEJ2000'] * u.deg
		positions_stars.append(SkyCoord(ra,dec,frame = 'icrs'))
		print dat['_2MASS']
		stars.append(dat['_2MASS'])

	return positions_stars,stars

#Open a LCO image fits file
def open_fits(filename):

	hdulist = fits.open(DIR+filename,ignore_missing_end=True)
	data = hdulist['SCI'].data
	header = hdulist['SCI'].header
	hdulist.close()

	return data,header

#Performs a cutout of the image around a single star
def cutout2D(data,position,header,x=50.0,y=50.0,unit = u.arcsec):
	size = u.Quantity((x,y), unit)
	cutout = Cutout2D(data, position, size, wcs=WCS(header))
	image_cutted = cutout.data

	return image_cutted


def find_stars(image_cutted,fwhm = 10.0):
	'''
	This function finds the centroid of the stars in the image above a threshold of 2 sigma from the noise level (median).
	:param: image_cutted: cutout of the full frame around the star of interest
	:param: fwhm: reference fwhm of the stars. This value doesn't need to be exact, can be approximated
	'''
	mean, median, std = sigma_clipped_stats(image_cutted,sigma = 3.0, iters = 20) 
	threshold = median + (2.0 * std)
	sigma_psf = float(fwhm)/gaussian_sigma_to_fwhm
	daofind = DAOStarFinder(fwhm=sigma_psf*gaussian_sigma_to_fwhm, threshold=threshold,sky = median) 
	#tbl = find_peaks(image_cutted, 2*threshold, box_size=5)

	sources = daofind(image_cutted-median)
	argsources = sources.argsort('flux') 
	argsources = np.flip(argsources,0)
	sources = sources[argsources]

	return sources


positions_stars, stars = choose_stars(stars_file,9)
print positions_stars

files = open_files(files_list)

#Loads the optimal aperture radii
rs = np.load(DIR_TO_SNR+'SNR_radius_dic_'+target_name+'.npy')

#This cycle performs aperture photometry on all files from all the selected stars. 
for star_ind in range(len(stars)):
	star_ind = 5
	appertures = []
	appertures_error = []
	date = []
	name = stars[star_ind]
	position = positions_stars[star_ind]
	means = []
	medians = []
	stds = []
	index = 0
	cameras = []
	telescopes = []
	seeings = []
	appertures_radius = []
	exptimes = []
	filters = []
	xs = []
	ys = []
	background = []
	radius = []
	file_index = 0
	centers_diago = []
	for filename in files:

		file_index += 1 
		print "File index = " + str(file_index)

		try:
			data,header = open_fits(filename)
			rdnoise = header['RDNOISE']
			mjd = header['MJD-OBS']
		except:
			print "Not working: " + str(filename)
			continue

		print "Working: " + filename
		print 'MJD: '  + str(mjd)
		print 'Star Name: ' + name
		image_cutted = cutout2D(data,position,header,50.0,50.0,u.arcsec)
		sources = find_stars(image_cutted,fwhm = header['L1FWHM'])

		positions = []

		r =  rs.item().get(filename) # Adaptative Aperture
		print "r = " +str(r)
		if r == None:
			print "No r measured"
			continue 

		#Some cameras have different fileds of view so we select the limits of where the centroid of the star of interest can be manually for each camera.   
		if True:
			if len(sources['xcentroid']) > 0:
				if True:
					for centroids in range(len(sources['xcentroid'])):
						if ('fl' in filename) and (sources['xcentroid'][centroids] < 69 ) & (sources['xcentroid'][centroids] > 60) & (sources['ycentroid'][centroids] > 60) & (sources['ycentroid'][centroids] < 69):
							source_x = sources['xcentroid'][centroids]
							source_y = sources['ycentroid'][centroids]
							positions.append([sources['xcentroid'][centroids], sources['ycentroid'][centroids]])
							print 'FL'
							break

						elif ('kb9' in filename) and (sources['xcentroid'][centroids] < 24 ) & (sources['xcentroid'][centroids] > 18) & (sources['ycentroid'][centroids] > 18) & (sources['ycentroid'][centroids] < 24):
							source_x = sources['xcentroid'][centroids]
							source_y = sources['ycentroid'][centroids]
							positions.append([sources['xcentroid'][centroids], sources['ycentroid'][centroids]])
							print 'KB9'
							break

						elif ('kb7' in filename) and (sources['xcentroid'][centroids] < 60 ) & (sources['xcentroid'][centroids] > 45) & (sources['ycentroid'][centroids] > 45) & (sources['ycentroid'][centroids] < 60):
							source_x = sources['xcentroid'][centroids]
							source_y = sources['ycentroid'][centroids]
							positions.append([sources['xcentroid'][centroids], sources['ycentroid'][centroids]])
							print 'KB7'
							break

						elif ('fs' in filename) and (sources['xcentroid'][centroids] < 88 ) & (sources['xcentroid'][centroids] > 80) & (sources['ycentroid'][centroids] > 80) & (sources['ycentroid'][centroids] < 88):
							source_x = sources['xcentroid'][centroids]
							source_y = sources['ycentroid'][centroids]
							positions.append([sources['xcentroid'][centroids], sources['ycentroid'][centroids]])
							print 'FS'
							break

					if len(positions) == 0: 
						#Uncomment to visually check the frames 
						#plot_image(image_cutted,header,str(exptime) + 'No star has been found')
						#source_x = sources['xcentroid'][0]
						#source_y = sources['ycentroid'][0]
						#positions.append([sources['xcentroid'][0], sources['ycentroid'][0]])
						print 'No star has been found'
						continue

				else:
					source_x = sources['xcentroid'][0]
					source_y = sources['ycentroid'][0]
					positions.append([sources['xcentroid'][0], sources['ycentroid'][0]])

				centers_diago.append(positions[0])

				if plot == True:
					fig, ax = plt.subplots()
					norm = ImageNormalize(image_cutted, interval=ZScaleInterval(),stretch=LinearStretch())
					plt.imshow(image_cutted, cmap='Greys', origin='lower', norm=norm)
					plt.plot(source_x,source_y,'+',color = 'red')
					cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(filename))
					plt.show()
					fig.canvas.mpl_disconnect(cid)
					plt.close(fig)

				index += 1
				print "We are in the sourcy file: " + str(index)

				r_in = 15
				r_out = 20 

				aper =  aper_phot(image_cutted,positions,r = r, r_in = r_in,r_out = r_out,bkg_sub = bkg_sub,plot=plot_aper)
				appertures.append(aper['residual_aperture_sum'][0])
				appertures_error.append(np.sqrt(aper['residual_aperture_sum'][0] + aper['aperture_sum_1'][0] + rdnoise**2))
				date.append(mjd)
				radius.append(r)
				background.append(aper['aperture_sum_1'][0])
				xs.append(source_x)
				ys.append(source_y)
				cameras.append(filename[-27:-23])
				telescopes.append(filename[0:4])
				seeings.append(header['L1FWHM'])
				filters.append(header['FILTER'])
				exptimes.append(header['EXPTIME'])
				#plot_aper = False
				medians.append(np.median(aper['aperture_sum_1']))
				means.append(np.mean(aper['aperture_sum_1']))
				stds.append(np.std(aper['aperture_sum_1']))

			else:
				print 'No stars found in the file: ' + filename	
		else:
			print 'No stars found in the file(except): ' + filename	

	r = '' #Adaptative Aperture
	print "Number of working files = " + str(len(appertures))

	#Centroid Distribution: Produces a 2D histogram with the distribution of the centroid. Usefull to see if the code is performing well. 
	centers_diago = np.array(centers_diago)
	centers_diago = centers_diago.astype(int)
	plt.hist2d(centers_diago[:,0],centers_diago[:,1],bins = 50,range = [[30,80],[30,80]], norm=colors.LogNorm())
	plt.colorbar()	
	plt.title(name + 'Centroid Distribution using DAOStarFinder')
	plt.savefig(DIR_TO_SAVE_FLUX+'DAO_Centroid_Distribution_mag.pdf')
	plt.clf()
	plt.show()

	#Plotting and saving the results
	if fluxes_study == True:
		appertures = np.array(appertures)
		appertures_error = np.array(appertures_error)
		date = np.array(date)
		exptimes = np.array(exptimes)
		cameras = np.array(cameras)
		telescopes = np.array(telescopes)
		seeings = np.array(seeings)
		xs = np.array(xs)
		ys = np.array(ys)
		radius = np.array(radius)
		background = np.array(background)

		plt.errorbar(date,appertures,yerr = appertures_error,fmt='o')
		plt.title('Flux ' + name)
		plt.xlabel('Date [MJD]')
		plt.ylabel('Flux [Counts]')
		plt.savefig(DIR_TO_SAVE_FLUX+'Flux_'+name+'_'+str(r)+'bkgsub.pdf')
		plt.close()
		np.save(DIR_TO_SAVE_FLUX+'Flux_'+name+'_'+str(r)+'bkgsub.npy',(date,appertures,appertures_error,xs,ys,cameras,seeings,exptimes,radius,background,filters))

	#Plotting and saving the results from the bakground study
	if bkg_study == True:
		plt.clf()
		plt.hist(np.array(medians))
		plt.title('Median ' + str(star_ind))
		plt.savefig(DIR_TO_SAVE_BACKGROUND + 'Median_' + str(name))
		plt.clf()
		plt.hist(np.array(means))
		plt.title('Mean ' + str(star_ind))
		plt.savefig(DIR_TO_SAVE_BACKGROUND + 'Mean_' + str(name))
		plt.clf()
		plt.hist(np.array(stds))
		plt.title('STD ' + str(star_ind))
		plt.savefig(DIR_TO_SAVE_BACKGROUND + 'STD_' + str(name))
		plt.clf()
		np.save(DIR_TO_SAVE_BACKGROUND + 'Anulus_Background' + str(name) + '.npy',(date,medians,means,stds))

	#Saving the file with the frames to skip 
	if bad_files == True: 
		files_not_use = np.array(files_not_use)
		np.save('bad_files_'+target_name+'.npy',files_not_use)