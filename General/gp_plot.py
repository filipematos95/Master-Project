import numpy as np
import matplotlib.pyplot as plt
import os,sys


"""
This code plots the results from the gp_fit.py. 
It outputs a plot with with lightcurve and the gp predictions within 1 sigma
"""

if len(sys.argv) > 2:
	period = sys.argv[1]
	if sys.argv[2] == 'B':
		DIR_TO_SAVE = #Write here the directory to save the output
		DIR_TO_LOAD = #Write here the folder with the data to load
		x,y,yerr = #Write here the file inside the DIR with the lightcurve
  	elif sys.argv[2] == 'O':
		DIR_TO_SAVE = #Write here the directory to save the output
		DIR_TO_LOAD = #Write here the folder with the data to load
		x,y,yerr = #Write here the file inside the DIR with the lightcurve
  	else:
  		print "You should choose a filter: B or R"
  		sys.exit()

x_pred,pred,pred_var = np.loadtxt(DIR_TO_LOAD+period+"_predictions.txt")
x,y,yerr,pred_r,pred_var_r = np.load(DIR_TO_LOAD+str(period)+"_data&pred.npy")

#Plot the fitted curve
plt.fill_between(x_pred, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var),color="k", alpha=0.2) #Plots 1 sigma error band
plt.plot(x_pred, pred, "k", lw=1.5, alpha=0.5)
plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)

if len(sys.argv) > 3:
	plt.xlim([np.min(x)-50,np.max(x)+50])
	plt.savefig(DIR_TO_SAVE+"GP_"+period[:6]+'_big.pdf')
else:
	plt.xlim([np.min(x)-20,np.max(x)+20])
plt.savefig(DIR_TO_SAVE+"GP_"+period[:6]+'.pdf')
plt.show()

