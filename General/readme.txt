Here are some instructions in order to use the codes, a more detailed documentation can be read in the my thesis: 

For all the codes: 
 - All the directories need to be set up. Please see the comments inside each script


get_wcs.py: 
 - It takes as input a text file with the names of the files to solve. 
 - It needs the client.py (they sould be in the same folder)
 - The command line to run should be: 
 	$ python get_wcs.py --newfits newname --wait --apikey [PERSONAL APIKEY]
 This way the output will be new fits files with updated headers.

SNR_study.py:
 - It needs the aperture.py module
Input: 
 - A fits file with the target star information. The header should contain: _2MASS,RAJ2000,DECJ2000,Bmag,Jmag 
 - The fits files with the LCO frames.
Output: 
 - A text file with the optimised aperture radius for each frame
 - A histogram with the distribution of the optimal redii
 - The light curve of the target star.

aperture_photometry.py: 
 - It needs the aperture.py module
 - Input:
	- A fits file with all the stars to be analised.The header should be : _2MASS,RAJ2000,DECJ2000,Bmag,Jmag
	- A text file with the name of all the frames to be analised.
        - All the frames mentioned above
 	- A text file with the optimised aperture radii computed by SNR_study.py 
 - Output: 
	- A light curve of all the stars before doing any differential photometry. 
	- 2D histogram with the distribution of the centroids (optional). 
	- A text file with all the files that should be skipped (optional).
 	- Pickle files with the measurements of the backgroud
	- A plot with the distribution of the background.

bkg_check.py:
 It takes as input the values of the background measured by aperture_photometry.py and produces 2D histograms of the background distribution of a full frame.

diferetential_photometry.py: 
 - Inputs: 
	- A fits file with the names of the stars. The header should be as the previous ones. 
	- The light curves produces by aperture_photometry.py (the directories should be set up)
- The user needs to select which reference stars wants to use.
 - Outputs: 
	- The final light curve of the target.
	- A plot with the RMS/photon noise of the produced light for every combination of reference stars.
	- A plot with the RMS for every combination of reference stars.
	- A text file with the legend for the combinations. Which stars are in which combination. 

stacking.py 
 - Stacks several images into one using the median. It produces a pdf with the image and a fits image. 
 - The command line should be $ python stacking.py [list] 
where list represents the name of a text file the names of the frames to stack. Alternatively the user can specify it inside the code. 

periodogram.py
 - It takes as input a time-series (x,y,yerr). 
 - It outputs a plot with the periodogram and a phase-folded light-curve (for the strongest peak) 

periodogram_bootstrap_multicore.py
 This code makes use of multicore machines. 
  - I takes as input a time-series (x,y,yerr), the number of cores to use and the number of periodograms to compute.
  -  It outputs a pickle file with a numpy array with the following: frequency, sum of LS power for all the periodograms, FAP,peaks
 
gp_fit.py
- Input:
	- The light curve (x,y,yerr) in which we want to train the model 
	- We can include two directories for light curves which can be selected in the command line. This are usefull when we have data from the same target in 2 filters or two seasons. 
	ex: $ python gp_fit.py B
- Outputs: 
	- initial hyperparameters guesses (text file) , optimized hyperparameters(text files) , mcmc distribution (pickle (npy) file) ,corner plot (pdf) and the gp predicitions(pickle (npy) and text file). 

gp_plot.py
This code is exclusive to plot the result from gp_fit.py
Input: 
 - The period of the computed model. $ python gp_plot.py xxx.xxxxxx filter
where xx.xxx represents the period in days (as outputed from gp_fit.py) and filter represents the B filter or R filter as introduced in the gp_fit.py