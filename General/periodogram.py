"""
This code produces a single periodogram of a light curve (or any time series)
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
rand = np.random.RandomState(42)
from astropy.stats import LombScargle
import sys
#plt.style.use('dark_background')

DIR = #Directory with the files to load
DIR_SAVE = #Directory to save the files

star_time, star_flux, star_flux_err =  np.load(DIR+'file_to_load')

print("There are: " + str(len(star_flux)))
print("Begin: " + str(np.min(star_time)))
print("End: " + str(np.max(star_time)))


ls = LombScargle(star_time, star_flux, star_flux_err)

if len(sys.argv) > 2: 
     min_period = sys.argv[1]
     max_period = sys.argv[2]
else:
     "You should specify the min and max periods. The default is being used (0.1:50 days)" 
     min_period = 0.1
     max_period = 50

freq, PLS = ls.autopower(minimum_frequency=1./max_period,
                         maximum_frequency=1./min_period)

best_freq = freq[np.argmax(PLS)]

phase = (star_time * best_freq) % 1

#Compute the best-fit model
phase_fit = np.linspace(0, 1)
mag_fit = ls.model(t=phase_fit / best_freq,
                   frequency=best_freq)
period = 1./best_freq

#Set up the figure & axes for plotting
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Lomb-Scargle Periodogram (period= ' + str(period)  + 'days)')
fig.subplots_adjust(bottom=0.12, left=0.07, right=0.95)
inset = fig.add_axes([0.78, 0.56, 0.15, 0.3])

#Plot the raw data
ax[0].errorbar(star_time, star_flux, star_flux_err, fmt='.', elinewidth=1.5, capsize=0)
ax[0].set(xlim=(np.min(star_time) -20, np.max(star_time)+20),
          xlabel='Observation time (days)',
          ylabel='Observed Flux')

#Plot the periodogram
ax[1].plot(1. / freq, PLS)
ax[1].set(xlabel='period (days)',
          ylabel='Lomb-Scargle Power',
          xlim=(0.1, 50.0),
          ylim=(0, 1));

#Plot the false-alarm levels
z_false = ls.false_alarm_level(0.001)
ax[1].axhline(z_false, linestyle='dotted', color='black')
print("Peak at: " + str(best_freq) + "\nFAP: " + str(ls.false_alarm_probability(PLS.max())))

#Plot the phased data & model in the inset
inset.errorbar(phase, star_flux, star_flux_err, fmt='.k', capsize=0)
inset.plot(phase_fit, mag_fit)
inset.invert_yaxis()
inset.set_xlabel('phase')
inset.set_ylabel('Normalized Flux [Counts]')
plt.savefig(DIR_SAVE+'filename to save')
plt.show()

#Plot the periodogram in a single figure
fig, ax = plt.subplots()
ax.plot(1. / freq, PLS)
ax.set(xlabel=r'Period (days)',
          ylabel=r'Lomb-Scargle Power',
          xlim=(0.1, 10.0),
          ylim=(0, 0.5))
ax.axhline(z_false, linestyle='dotted',color = 'black')
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')
ax.xaxis.label.set_color('black')
ax.yaxis.label.set_color('black')

ax.grid(color='grey', linestyle='--', linewidth=0.5)

#ax.set_facecolor('xkcd:white')
plt.savefig(DIR_SAVE+"filename to save",dpi=500)
plt.show()
