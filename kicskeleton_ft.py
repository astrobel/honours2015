import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import smoothing

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{helvet}',
       r'\usepackage[EULERGREEK]{sansmath}',
       r'\sansmath'
]

def fft(time, flux, ofac, hifac):
  # Do LNP Test (Lomb-Scargle Periodogram) ie FT
  freq,power, nout, jmax, prob = lomb.fasper(time, flux, ofac, hifac)
  convfactor = (1. / (60 * 60 * 24)) * (10 ** 6)
  uHzfreq = freq * convfactor #11.57, conversion c/d to mHz

  return uHzfreq, power


# read in light curve
importblend = np.loadtxt('kicskeleton_lc.dat')
clipped_time = importblend[:,0]
clipped_flux = importblend[:,1]

# fourier transform
print ' '
print '~~~Fourier transform for power spectrum~~~'
ofac = raw_input('Oversampling factor (low): ')
hifac = raw_input('Hifac: ')
print ' '
ofac = np.float64(ofac)
hifac = np.float64(hifac)
frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
hifac = 283 / max(frequencies)
frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
if power_spectrum.size != frequencies.size:
   power_spectrum = np.concatenate((power_spectrum, [0]))
power_spectrum = power_spectrum * 4 * np.var(clipped_flux) / clipped_flux.size
power_spectrum = np.sqrt(power_spectrum)
power_spectrum *= 1e6 #* ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

# logarithmic and linear power spectra
fig1, (pslog, pslin) = plt.subplots(2, 1) 
pslog.plot(frequencies, power_spectrum, 'r-')
pslog.set_xlim(1, max(frequencies))
pslog.set_ylabel('Amplitude (ppm)')
pslog.set_xscale('log')
pslog.set_yscale('log')
pslog.set_title('skeleton')

pslin.plot(frequencies, power_spectrum, 'r-')
pslin.set_xlim(1, max(frequencies))
pslin.set_ylim(ymin = 0)
pslin.set_xlabel('Frequency ($\mu$Hz)')
pslin.set_ylabel('Amplitude (ppm)')

plt.tight_layout()
fig1.savefig('kicskeleton_spectrum.png')

# plt.show()
