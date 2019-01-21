'''
This code compares the different apperture photometry appertures and stars by comparing the RMS of the errors.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from itertools import combinations
import warnings
from astropy.io import fits
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

DIR = #Directory with the files to load
DIR_TO_SAVE = #Directory to save the files
plot = False #True to show interactive plots of each image while running

def star_fits (filename,sort):
	"""
	Function to open fits files with the list of stars
	:params: filename: Filename of the fits file with the list of stars
	:params: sort: Collumn name of the collumn to sort (ex: Bmag)
	:out: data: Array witht the sorted stars information 
	"""

	hdulist = fits.open(filename,ignore_missing_end=True)
	data = hdulist[1].data
	data = data[data[sort].argsort()]

	return data

def outlier_removal(dataset,k,ite=1):

	"""
	Function to removes the outliers. It removes data points that are k rms from the median
	:params: dataset: time series (y values)
	:params:ite: Number of iterarions
	:out: dataset_clean: time seris (y values) without the outliers
	:out: len(removed[0]): Number of outliers removed
	"""

	removed = [[1]]
	for i in range(ite):
		median = np.median(dataset)
		residuals = dataset - median
		rms = np.sqrt(np.mean(np.square(residuals)))
		dataset_clean = np.where(np.absolute(residuals) <= (k*rms))
		removed = np.where(np.abs(residuals) > k*rms)
		dataset = dataset[dataset_clean]

	return dataset_clean,len(removed[0])


datas = star_fits('filename with the stars list.fit','Bmag')
positions_stars = []
stars_list = []

#Load the stars to analyze
for dat in datas[0:10]: #The number of stars should be chosen here
	#if dat['_2MASS'] != 'NAME OF THE SATAR TO EXCLUDE': #Exclude the target star for example
	stars_list.append(dat['_2MASS'])

stars_list = np.array(stars_list)
fluxes = np.array([]) 
fluxes_day = np.array([])
r = 11 #Default value of the aperture radius (only used in case that there is no particular optimisation based on the SNR)
rms_list = []
indo = 0

for star1 in stars_list:
	text_file = open(DIR_TO_SAVE+'Combiations_'+str(star1)+'_'+str(r)+'.txt', 'w')
	combination = []
	star1_time, star1_flux,star1_flux_err = np.load(DIR+'Flux_'+star1+'_'+'_out_error.npy')

	#Normalizes de flux
	star1_flux_norm = star1_flux/(np.median(star1_flux))
	star1_flux_err = (star1_flux_err/star1_flux)*star1_flux_norm
	star1_flux = star1_flux_norm

	print 'Flux ' + str(star1) + ' = ' + str(np.sum(star1_flux)) 
	stars_list1 = stars_list[np.where(stars_list != star1)]
	rms_list = []
	indo = 0
	for n_stars in range(1,len(stars_list1)+1):
		for star_list in list(combinations(stars_list1,n_stars)):
			star_time = star1_time
			star_flux = np.zeros(len(star1_time))
			star_flux_err = np.zeros(len(star1_time))

			for star2 in star_list:
				star2_time, star2_flux ,star2_flux_err = np.load(DIR+'Flux_'+star2+'_'+'_out_error.npy')
				
				#Normalizes the flux
				star2_flux_norm = star2_flux/(np.median(star2_flux)) 
				star2_flux_err = (star2_flux_err/star2_flux)*star2_flux_norm
				
				star2_flux = star2_flux_norm[np.in1d(star2_time,star_time)]

				star2_flux_err = star2_flux_err[np.in1d(star2_time,star_time)]
				star_flux_err = star_flux_err[np.in1d(star_time,star2_time)]
				star_flux = star_flux[np.in1d(star_time,star2_time)]
				star_time = star_time[np.in1d(star_time,star2_time)]
				star_flux += star2_flux

				star_flux_err += (star2_flux_err/star2_flux)**2
			
			#Normalizes de flux
			star_flux_norm = star_flux/(np.median(star_flux))
			star_flux_err = (star_flux_err/star_flux)*star_flux_norm
			star_flux = star_flux_norm

			star1_flux = star1_flux[np.in1d(star1_time,star_time)]
			star1_flux_err = star1_flux_err[np.in1d(star1_time,star_time)]
			star1_time = star1_time[np.in1d(star1_time,star_time)]
			rel_photometry = star1_flux/star_flux

			erro = rel_photometry*np.sqrt( (star1_flux_err/star1_flux)**2 + star_flux_err) 
			rel_photometry_error = erro
			rel_photometry_error = rel_photometry_error[~np.isnan(rel_photometry_error)]
			rel_photometry = rel_photometry[~np.isnan(rel_photometry_error)]
			star_time = star1_time[~np.isnan(rel_photometry_error)]

			removed = 0 #0 to disable outliers removal 

			#Remove outliers 
			while removed > 0:
				no_outliers,removed = outlier_removal(rel_photometry,3)

				rel_photometry = rel_photometry[no_outliers]
				star_time = star_time[no_outliers]
				rel_photometry_error = rel_photometry_error[no_outliers]


			median = np.median(rel_photometry)
			residuals = rel_photometry - median
			
			print 'Combination: ' + str(indo) #Printing status 
			print star_list		
			text_file.write(str(indo)+' '+str(star_list)+'\n')
			combination.append(stars_list)

			rms = np.sqrt(np.mean(np.square(residuals)))
			mean_error = np.mean(rel_photometry_error)

			plt.errorbar(star_time,rel_photometry,yerr = rel_photometry_error,fmt='.')
			plt.axhline(y=median,color='red',linestyle='--')	
			plt.errorbar(((star_time[-1] + star_time[0] )/ 2.),median,yerr = rms, fmt = '.')
			plt.title('Residuals: ' +  str(indo) + ' ' + str(star1))
			#plt.ylim([-4*rms,4*rms])
			r = 'r_best'
			plt.savefig(DIR_TO_SAVE+str(r)+'/Residuals_'+str(r)+'_'+str(indo)+'_'+str(star1)+'.pdf') 
			#plt.show()
			plt.clf()
			np.save(DIR_TO_SAVE+str(r)+'/Residuals_'+str(r)+'_'+str(indo)+'_'+str(star1)+'.npy',(star_time,rel_photometry,rel_photometry_error))
			rms_list.append([star_list,rms,mean_error])
			indo += 1
	

	r = 'r_best'
	text_file.close()
	combination = np.array(combination)

	#Save the results
	np.save(DIR_TO_SAVE+str(r)+'/Combiations_'+str(star1)+'_'+str(r)+'.npy',combination)
	
	#Plotting
	rms_list = np.array(rms_list)
	rms_indices = range(len(rms_list))

	plt.plot(rms_indices,rms_list[:,1])
	plt.xlabel('Combination')
	plt.ylabel('Root mean square from Median')
	plt.title('Apperture: ' + str(r) + ' ' + str(star1))
	plt.savefig(DIR_TO_SAVE+str(r)+'/RMS_10'+ str(r)+'_'+ str(star1) +'_.pdf')
	plt.show()
	
	plt.clf()
	plt.plot(rms_indices,rms_list[:,1]/rms_list[:,2])
	plt.xlabel('Combination')
	plt.ylabel('Root mean square from Median / Photon Noise')
	plt.title('Apperture: ' + str(r) +' '+ str(star1))
	plt.savefig(DIR_TO_SAVE+str(r)+'/RMS'+ str(r)+'_'+str(star1)+'.pdf')
	plt.show()
	np.save(DIR_TO_SAVE+'NORM/'+str(r)+'/RMS'+ str(r)+'_'+str(star1)+'.pdf',rms_list)