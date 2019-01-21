"""
This code stacks several frames into one image. It computes the median for each pixel.
"""
from astropy.io import fits
import numpy as np 
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize,LinearStretch,ZScaleInterval)


DIR = #Directory with the frames to stack 

if len(sys.argv) < 2: 
	source = 'list.txt' #list.txt should contain the name of the frames to stack
else:
	source = sys.argv[1]

with open(DIR+source, 'r') as fp:
    files = fp.read().splitlines()

final_image = []
for filename in files[0:50]:
	hdulist = fits.open(DIR+filename,ignore_missing_end=True)
	image = hdulist['SCI'].data

	final_image.append(image)
	hdulist.close()

final_image = np.median(final_image,0)
data = final_image
norm = ImageNormalize(data, interval=ZScaleInterval(),stretch=LinearStretch())
plt.imshow(data, cmap='Greys', origin='lower', norm=norm)

hdr = fits.open(DIR+files[33],ignore_missing_end=True)['SCI'].header

fits.writeto('stacked_0_50.fits',final_image,hdr, overwrite=True)
plt.savefig('stacked_0_50.pdf',overwrite=True)
plt.show()