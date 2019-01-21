import george 
from george.modeling import Model
from george import kernels
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import os,sys

''' 

This code models the light curve using gaussian processes. 
The hyperparameter optimzation is done by sampling & maginalization using MCMC (emcee)

'''

#Importing the light curve data: The directory should be adapted for each user
#DIR_TO_SAVE: Where the output files will be saved 
#DIR: Where the light-curve/timeseries is located

if len(sys.argv) > 1:
	if sys.argv[1] == 'B':
		DIR_TO_SAVE = #Write here the directory to save the output
		DIR = #Write here the folder with the data to load
		x,y,yerr = #Write here the file inside the DIR with the lightcurve
  	elif sys.argv[1] == 'R':
		DIR_TO_SAVE = #Write here the directory to save the output
		DIR = #Write here the folder with the data to load
		x,y,yerr = #Write here the file inside the DIR with the lightcurve
  	else:
  		print "You should choose a filter: B or R"
  		sys.exit()

if not os.path.exists(DIR_TO_SAVE):
	os.makedirs(DIR_TO_SAVE) 

#Normalize de flux if its not normalized (optional)
if np.median(y) != 0:
	y = y - np.median(y)


class SinModel(Model):
	'''
	This class implements a mean model for the GP object with the type of a Sin wave
	:params: Model: GP object
	'''
	parameter_names = ("amp","per")

	def get_value(self,t):
		t = t.flatten()
		return (self.amp*np.sin(2*np.pi/self.per*t))


''' Mean Model '''

truth = dict(amp = 0.03, per = 4.43) #Mean Initial Parameters 
kwargs = dict(**truth)
kwargs["bounds"] = dict(amp = (0.03,0.08), per = (3.0,7.0)) #Mean Hyperparameters Bounds
mean_model = SinModel(**kwargs)

''' Kernels '''

bounds = dict(gamma = (0.01,25.0),log_period = (np.log(3.4),np.log(7.0))) #Bounds for the ExpSine2Kerner
kernel_gprot= np.var(y)*kernels.ExpSquaredKernel(5.0)*kernels.ExpSine2Kernel(gamma=5.1, log_period = np.log(4.43),bounds=bounds) #Kernel Definition
kernel_labels = ["amp","period","var","gamma_exp","gamma_sin","log_period"]
model = george.GP(kernel_gprot,fit_white_noise = True,mean=mean_model,fit_mean = True)

model.compute(x,yerr) #Initializes de model

print("Inital parameters:")
print(model.get_parameter_vector())
print(model.get_parameter_bounds())

inital_params = model.get_parameter_vector()
initial_bounds = model.get_parameter_bounds()

print(len(kernel_labels))
print(len(inital_params))
print(len(initial_bounds))
#Save the inital parameters
save_file_init = open(DIR_TO_SAVE+str(65)+"_init_hyperparams.txt", "w")
for i in range(len(kernel_labels)):
	save_file_init.write(str(kernel_labels[i])+": "+str(inital_params[i]) + ' ' + str(initial_bounds[i]) + '\n')
save_file_init.close()

#Weight function 
def lnprob(p):
    model.set_parameter_vector(p)
    return model.log_likelihood(y, quiet=True) + model.log_prior()

#MCMC for parameter optimization
initial = model.get_parameter_vector()
ndim, nwalkers = len(initial), 32
p0 = initial + 1e-6 * np.random.randn(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
print("Initial ln-likelihood: {0:.2f}".format(model.log_likelihood(y)))
print("Running burn-in...")
p0, _, _ = sampler.run_mcmc(p0, 500) #Define the number of samples to take in the burn-in
#sampler.reset()

print("Running production...")
sampler.run_mcmc(p0,1000); #Define the number of samples

samples = sampler.flatchain
model.set_parameter_vector(np.percentile(samples,[50], axis=0)[0])

#Compute the predictions
model.recompute()
print("\nFinal ln-likelihood: {0:.2f}".format(model.log_likelihood(y)))
period = model.get_parameter_vector()[1]
print("Final period: " +str(period))

n_days = np.max(x)-np.min(x)
x_pred = np.linspace(np.min(x)-50,np.max(x)+50,n_days*100)
pred, pred_var = model.predict(y, x_pred, return_var=True)
pred_r, pred_var_r = model.predict(y, x, return_var=True)

r_xi = np.sum((y-pred_r)**2/yerr**2) / (len(y)-4)
print "Chi: " + str(r_xi)

#Save predictions
np.savetxt(DIR_TO_SAVE+str(period)+"_predictions.txt",(x_pred,pred,pred_var))
np.save(DIR_TO_SAVE+str(period)+"_data&pred",(x,y,yerr,pred_r,pred_var_r))
np.save(DIR_TO_SAVE+str(period)+"_mcmc",(samples))

''' CORNER ''' #Produces the Corner Plot
print model.get_parameter_vector()
corner.corner(samples,labels = kernel_labels,quantiles=[0.16, 0.5, 0.84],show_titles=True, title_kwargs={"fontsize": 12});

plt.savefig(DIR_TO_SAVE+'GP_'+str(period)+'_corner.pdf')
plt.clf()


#Save the inital parameters
save_file_init = open(DIR_TO_SAVE+str(period)+"_init_hyperparams.txt", "w")
for i in range(len(kernel_labels)):
	save_file_init.write(str(kernel_labels[i])+": "+str(inital_params[i]) + ' ' + str(initial_bounds[i]) + '\n')
save_file_init.close()

mean = np.percentile(samples,[50], axis=0)[0]
plus = np.percentile(samples,[84], axis=0)[0]
minus = np.percentile(samples,[16], axis=0)[0]

#Save the optimized hyperparameters
save_file = open(DIR_TO_SAVE+str(period)+"_opt_hyperparams.txt", "w")
for i in range(len(kernel_labels)):
	save_file.write(str(kernel_labels[i])+": "+str(mean[i]) + '	+' + str(plus[i]) + ' -' + str(minus[i]) +'\n')
save_file.write("Chi: " + str(r_xi))
save_file.close()


