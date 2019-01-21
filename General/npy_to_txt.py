import numpy as np 
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
	file_name = sys.argv[1]

x_pred,pred,pred_var = np.loadtxt(file_name+"_predictions.txt")
x,y,yerr,pred_r,pred_var_r = np.load(file_name+"_data&pred.npy")
samples = np.load(file_name+"_mcmc.npy")

print(samples.shape)
#Plot the fitted curve
plt.fill_between(x_pred, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var),color="k", alpha=0.2) #Plots 1 sigma error band
plt.plot(x_pred, pred, "k", lw=1.5, alpha=0.5)
plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)

if len(sys.argv) > 3:
	plt.xlim([np.min(x)-50,np.max(x)+50])
else:
	plt.xlim([np.min(x)-20,np.max(x)+20])

np.savetxt(file_name+"_predictions_ascii.txt",np.transpose([x_pred,pred,np.sqrt(pred_var)]),header = 'Time [MJD],Prediction (Relative Flux), Error')
np.savetxt(file_name+"_ligth_curve_ascii.txt",np.transpose([x,y,yerr]),header = 'Time [MJD], Relative FLux, Error')
np.savetxt(file_name+"_mcmc_ascii.txt",(samples))