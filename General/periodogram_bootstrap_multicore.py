""""
This code estimates the errors of Lomb-Scargle Periodograms using a bootstrap approach.
This code allows the user to use multicores to minimize the computation time. There is no limit in the number of codes to use. 
"""

import multiprocessing as mp
import logging
import traceback
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('ggplot') 
from astropy.stats import LombScargle
import os
import sys
from scipy.signal import argrelextrema
from tqdm import tqdm
rand = np.random.RandomState(42)

info = mp.get_logger().info

def plot(freq,PLS,ls):
"""
Plotting function

:params: freq: array with the frequencies
:params: PLS: array with the LS powers
:params: ls: Lomb-Scargle Astropy obejct
"""

  fig, ax = plt.subplots()
  ax.plot(1. / freq, PLS,color='black')
  z_false = ls.false_alarm_level(1.0)
  ax.axhline(z_false, linestyle='dotted', color='black')
  ax.set(xlabel='period (days)',
          ylabel='Lomb-Scargle Power',
          xlim=(2, 20.0),
          ylim=(0, 0.3));
  ax.grid(color='grey', linestyle='--', linewidth=0.5)
  ax.set_facecolor('xkcd:white')
  plt.show()

def main(): 
    """
    Main function
    """

    DIR = #Direcotry to load the input light curve
    DIR_SAVE = #Directory to save the out files 
    pip = 'b' #Season of the data (based on the pipeline b (BANZAI - 2nd Season) or o (ORAC - 1st Season)) 
    star_time, star_flux, star_flux_err =  np.load(DIR+'Residuals_r_best_8_55555555-0208391.npy')

    time = star_time
    flux = star_flux
    error = star_flux_err 

    if len(sys.argv) > 2:
        min_freq = 1/float(sys.argv[2])
        max_freq = 1/float(sys.argv[1])
        min_per = sys.argv[1]
        max_per = sys.argv[2]
    else:
        print("Please provide the limits as arguments [min] [max] \n Default is being used [0:10]")
        min_freq = 1/10.0
        max_freq = 1/0.01
        min_per = 1.0
        max_per = 10.0

    nproc = mp.cpu_count()  
    nproc = max(1, nproc) - 1

    logger = mp.log_to_stderr()
    logger.setLevel(logging.INFO)

    ntasks = 3 #Number of cores to use
    n_samples = 33000 #Number of samples to use 
    inputs = [[time,flux,error,n_samples,min_freq,max_freq],[time,flux,error,n_samples,min_freq,max_freq],[time,flux,error,n_samples,min_freq,max_freq]]

    #inputs = [(,) for i in xrange(ntasks)]  #Tuple with inputs for each worker

    in_q  = mp.Queue()
    out_q = mp.Queue()
    procs = [mp.Process(target=worker, args=(in_q, out_q)) for i in range(nproc)]

    for i in range(ntasks):
        in_q.put(inputs[i])

    for i in range(nproc):
        in_q.put('STOP')

    for p in procs:
        p.start()

    allout = []
    peaks = np.array([])
    per_sum = 0

    while ntasks > 0:
        result = out_q.get() #Get data from workers
        peaks = np.concatenate((peaks,result[1]))
        if type(per_sum) is int:
            per_sum = result[0]
        else: 
            per_sum += result[0]

        ntasks -= 1

    for p in procs:
        p.join()

    ls = LombScargle(time, flux, error)
    fal = ls.false_alarm_level(0.5)
    freq, PLS = ls.autopower(minimum_frequency=min_freq,
                        maximum_frequency=max_freq,samples_per_peak=1000)
    plot(freq,per_sum/(n_samples*3),ls)
    plt.hist(1/peaks,200)
    plt.show()
    np.save(DIR_SAVE+"bootstrap_"+str(min_per)+"_"+str(max_per)+"_"+str(pip)+"_samp_"+str(n_samples*3)+".npy",(freq,per_sum,fal,peaks))
    return 0

def worker(in_q, out_q):
    """
    Work done in each child process.
    :params: in_q: is the input Queue
    :params: out_q is the output Queue.
    """
    while True:
        try:
            tmp = in_q.get()
            if tmp == 'STOP':
                break
            out_dict = {}

            out_dict = bootstrap(*tmp) #Do stuff here

            out_q.put(out_dict)

        except Exception as exception:
            info(str(traceback.format_exc()))
            return
    return 

def bootstrap(time,flux,error,size,min_freq,max_freq):
    """
    Function that performs the bootstraping
    :params: time: time series array (x value)
    :params: flux: time series array (y value)
    :params: error: error on the y 
    :parms: size: 
    :params: min_freq: minimum frequency to search 
    :params: min_freq: maximum frequency to search

    :out: per_sum: array with the sum of all the Periodograms
    :out: peaks: array with frequency of each peak 
    """


    flux_distributions = np.zeros((len(flux),size))
    for i in range(0,len(flux)):
        flux_distributions[i] = np.random.normal(flux[i],error[i],size)

    peaks = []
    per_sum = 0  
    peaks_len = 0
    for i in range(size):
        ls = LombScargle(time, flux_distributions[:,i], error)
        threshold = ls.false_alarm_level(0.99)
        freq, PLS = ls.autopower(minimum_frequency=min_freq,maximum_frequency=max_freq,samples_per_peak=1000)
        if type(per_sum) is int:
            per_sum = PLS
        else: 
            per_sum += PLS

        peaks.append(freq[argrelextrema(PLS,np.greater)])

        if i%(size/100) == 0:
            print("Progress: %i %%"%(i/size*100))

    peaks = np.concatenate(peaks)
    return(per_sum,peaks)

if __name__ == "__main__":
    main()
