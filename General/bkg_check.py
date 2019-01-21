"""
This script checks the background distribution of fits images
"""

import numpy as np
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils import Background2D, MedianBackground
from photutils.background import MMMBackground, MADStdBackgroundRMS
import time
DIR = #Directory with the files to load
DIR_TO_SAVE = #Directory to save the output files

with open(DIR+'list.txt', 'r') as fp: #Opens the list with files to analise. The file should have the name list.txt
	files = fp.read().splitlines()

	print len(files)
	
	means = []
	medians = []
	stds = []
	means2 = []
	medians2 = []
	stds2 = []
	means = np.zeros(shape=(3,3))
	medians = np.zeros(shape=(3,3))
	stds = np.zeros(shape=(3,3))
	size_max = 3999
	size = 1333

	for filename in files[200:-1]: 

		try:
			hdulist = fits.open(DIR+filename,ignore_missing_end=True)
		except:
			print filename
			continue
		print "Working: " + filename

		data = hdulist[0].data
		data = data[11:4010,11:4010]

		#mean_, median_, std_ = sigma_clipped_stats(data,sigma = 3.0, iters = 5)
	
		blocks = map(lambda x : np.split(x, data.shape[1]/size, 1),np.split(data, data.shape[0]/size, 0))
		
		for i in range(size_max/size):
			for j in range(size_max/size):
				mean, median, std = sigma_clipped_stats(blocks[i][j],sigma = 3.0, iters = 5)
				means[i][j] += mean
				medians[i][j] += median
				stds[i][j] += std
				means2.append((mean-mean_)/mean_)
				medians2.append((median-median_)/median_)
				stds2.append((std-std_)/std_)

	
	#Plots the distribution of tha noise in the image 

	plt.hist(means2)
	plt.grid(True)
	plt.show()
	plt.hist(medians2)
	plt.grid(True)
	plt.show()
	plt.hist(stds2)
	plt.grid(True)
	plt.show()
	
	np.save(DIR_TO_SAVE+'bkg_check_total33_3.npy',(means,medians,stds))


	#Plots 2D histrograms with the spacial distribution around the images. 
	plt.imshow(means)
	plt.colorbar()
	plt.show()
	plt.imshow(medians)
	plt.colorbar()
	plt.show()
	plt.imshow(stds)
	plt.colorbar()
	plt.show()


