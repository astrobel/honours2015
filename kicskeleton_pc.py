import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import smoothing
import translate as tr

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{helvet}',
       r'\usepackage[EULERGREEK]{sansmath}',
       r'\sansmath'
]
mpl.rcParams['axes.formatter.useoffset'] = False

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

# second, more intense fourier transform
print ' '
print '~~~Fourier transform for phase curve~~~'
ofac = raw_input('Oversampling factor (high): ')
hifac = raw_input('Hifac: ')
topfreq = raw_input('Highest frequency to plot in microHertz: ')
print ' '
ofac = np.float64(ofac)
hifac = np.float64(hifac)
topfreq = np.float64(topfreq)
frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
hifac = topfreq / max(frequencies)
frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
power_spectrum = np.concatenate((power_spectrum, [0]))
power_spectrum = power_spectrum * 4 * np.var(clipped_flux) / clipped_flux.size
power_spectrum = np.sqrt(power_spectrum)
power_spectrum *= 1e6 #* ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

# close-up of interesting part of power spectrum
plt.figure(1)
plt.plot(frequencies, power_spectrum, 'r-')
plt.xlim(0, max(frequencies))
plt.ylim(ymin = 0)
plt.xlabel('Frequency ($\mu$Hz)')
plt.ylabel('Amplitude (ppm)')
plt.title('skeleton')
plt.savefig('kicskeleton_zoomed.png')

# detecting frequency for phase curve and folding
foldfreq = frequencies[power_spectrum.argmax()]
convfactor = (1. / (60 * 60 * 24)) * (10 ** 6)
foldfreq = foldfreq / convfactor
foldper = 1. / foldfreq
clipped_time = clipped_time % foldper

# sort time and flux for binning
sortblend = np.zeros([clipped_time.size, 2])
sortblend = np.array([clipped_time, clipped_flux])
np.reshape(sortblend, (2, clipped_time.size))
sortblend.sort(axis=0)
np.reshape(sortblend, (clipped_time.size, 2))
sortblend = np.transpose(sortblend)
time = sortblend[:,0]
sap_flux = sortblend[:,1]

smoothtype2 = raw_input('Smooth or bin the phase curve? (s/b) ')
while smoothtype2 != "b" and smoothtype2 != "s":
   smoothtype2 = raw_input('Please type either b or s: ')
print ' '

if smoothtype2 == "s": # smoothing

   print '~~~Phase curve smoothing~~~'
   smoothtype = raw_input('Boxcar or Gaussian smoothing? (b/g) ')
   while smoothtype != "b" and smoothtype != "g":
      smoothtype = raw_input('Please type either b or g: ')
   kern = raw_input('Smoothing kernel: (integer) ')
   kern = np.float64(kern)
   print ' '

   # smoothing method
   if smoothtype == "b": # boxcar smoothing
      sap_flux2, smth_flux = smoothing.boxsmooth(time, sap_flux, kern)
   elif smoothtype == "g": # gaussian smoothing
      sap_flux2, smth_flux = smoothing.gausssmooth(time, sap_flux, kern)

   # plotting phase curve
   plt.figure(2)
   plt.plot(time, smth_flux, 'ro', markersize=3)
   plt.xlabel('Time mod %f days' % foldper)
   plt.ylabel('Fractional Intensity')
   plt.title('skeleton')
   plt.savefig('kicskeleton_phase.png')

elif smoothtype2 == "b": # binning

   print '~~~Phase curve binning~~~'
   binnum = raw_input('Number of bins: (integer) ')
   binnum = np.float64(binnum)
   print ' '

   # manual histogram
   binsize = max(time) / binnum
   bindex = 0
   flux_sum = np.zeros(max(time) / binsize)
   flux_num = np.zeros(max(time) / binsize)
   for i, t in enumerate(time):
      bindex = (t - (t % binsize)) / binsize - 1
      flux_sum[bindex] = sap_flux[i] + flux_sum[bindex]
      flux_num[bindex] += 1
   flux_sum = np.divide(flux_sum, flux_num)
   time_binned = np.linspace(0, max(time), binnum)

   time_doubled = np.zeros(time_binned.size)

   for i, val in enumerate(time_binned):
      time_binned[i] = tr.translate(val, 0, foldper, 0, 1)
      time_doubled[i] = time_binned[i] + 1

   binnum2 = binnum / 10.
   binsize2 = max(time) / binnum2
   bindex2 = 0
   flux_sum2 = np.zeros(max(time) / binsize2)
   flux_num2 = np.zeros(max(time) / binsize2)
   for i, t2 in enumerate(time):
      bindex2 = (t2 - (t2 % binsize2)) / binsize2 - 1
      flux_sum2[bindex2] = sap_flux[i] + flux_sum2[bindex2]
      flux_num2[bindex2] += 1
   time_binned2 = np.linspace(0, max(time), binnum2)
   flux_sum2 = np.divide(flux_sum2, flux_num2)

   time_doubled2 = np.zeros(time_binned2.size)

   for i, val in enumerate(time_binned2):
      time_binned2[i] = tr.translate(val, 0, foldper, 0, 1)
      time_doubled2[i] = time_binned2[i] + 1

   # plotting phase curve
   plt.figure(2)
   plt.plot(time_binned, flux_sum, 'ro', markersize=3)
   plt.plot(time_doubled, flux_sum, 'ro', markersize=3)
   plt.plot(time_binned2, flux_sum2, 'cs')
   plt.plot(time_doubled2, flux_sum2, 'cs')
   plt.xlim(0, max(time_doubled))
   plt.xlabel('Normalised Time mod %f days' % foldper)
   plt.ylabel('Fractional Intensity')
   plt.title('skeleton')
   plt.savefig('kicskeleton_phase.png')

# plt.show()
