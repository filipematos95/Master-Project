"""
Script to load the header information from LCO fits files
"""

from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
import sys

DIR = #Directory with the files to load

if len(sys.argv) < 2:
	print 'There is no file to open'
	sys.exit()

for i in range(1,len(sys.argv)):
	filename = DIR+sys.argv[i]
	print filename[-5:]
	if filename[-5:] != '.fits':
		print 'The file introduced is not a fits file'
	else:
		hd = fits.open(filename)
		print hd.info()

		header_number = raw_input("Which header do you wanna see?")
		header = hd[int(header_number)].header
		data = hd[int(header_number)].data
		#print repr(header)
		header_list = list(header.keys())  
		item = raw_input("Which item do you wanna see?")
		i_right = 0
		for i in range(0,len(header_list)):
			if header[i] == item:
				i_right = i
				break
		print header[i_right]



