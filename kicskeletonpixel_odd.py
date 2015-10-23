import numpy as np
import scipy as sp
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import lomb
import smoothing 
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl

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


### INPUT PARAMETERS ###

print ' '
print '~~~Light curve cleaning~~~'
smoothtype = raw_input('Boxcar or Gaussian smoothing? (b/g) ')
if smoothtype != "b" and smoothtype != "g":
   smoothtype = raw_input('Please type either b or g: ')
kern = raw_input('Smoothing kernel: (integer) ')
kern = np.float64(kern)

print ' '
inp = raw_input('Clipping level: ')
print ' '
inp = np.float64(inp)

print ' '
print '~~~Fourier transform~~~'
ofac = raw_input('Oversampling factor (low): ')
hifac = raw_input('Hifac: ')
print ' '
ofac = np.float64(ofac)
hifac = np.float64(hifac)


### REPRESENTATIVE QUARTER: 13 ###

hdulist = pyfits.open('kplr00skeleton-2012179063303_lpd-targ.fits', mode='readonly')

table = hdulist[1].data
flux = table['FLUX']
#flux_err = table['FLUX_ERR']
time = table['TIME']
hd1 = hdulist[1].header
ysize = hd1['NAXIS2']

table2 = hdulist[2].data
hd2 = hdulist[2].header
x = hd2['NAXIS1']
y = hd2['NAXIS2']
xsize = x * y
temp2d = np.zeros((x, y))

hdulist.close()

# dynamic variable names
for (j, k), img in np.ndenumerate(temp2d):
   index = (k + 1) * j + (x - j) * k
   exec("pixel%d_flux = np.array(None)" % index)
   exec("pixel%d_time = np.array(None)" % index)

second_flux = np.zeros([xsize, ysize])

# filling the flux array
for (i, j, k), val in np.ndenumerate(flux):
   index = (j + 1) * k + (x - k) * j
   second_flux[index, i] = val

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      pass
   else:
      flux2 = second_flux[index,:]

      # creating a blend array to remove NaNs
      blend = np.array([time, flux2])
      blend = np.transpose(blend)
      blend2 = np.ma.compress_rows(np.ma.fix_invalid(blend))

      time2 = blend2[:,0]
      flux2 = blend2[:,1]

      if smoothtype == "b": # boxcar smoothing
         flux3, smth_flux = smoothing.boxsmooth(time2, flux2, kern)
      elif smoothtype == "g": # gaussian smoothing
         flux3, smth_flux = smoothing.gausssmooth(time2, flux2, kern)

      exec("pixel%d_flux = flux3" % index)
      exec("pixel%d_time = time2" % index)

      exec("tempflux = pixel%d_flux" % index)
      exec("temptime = pixel%d_time" % index)

      clip = inp * np.std(tempflux)
      meanflux = np.mean(tempflux)
 
      upperbound = meanflux + clip
      lowerbound = meanflux - clip

      colours = np.zeros(tempflux.size)

      for i, flux in enumerate(tempflux):
         if flux < upperbound and flux > lowerbound:
            colours[i] = 1

      clipped_flux = []
      clipped_time = []
      for i, colour in enumerate(colours):
         if colour == 1:
            clipped_flux.append(tempflux[i])
            clipped_time.append(temptime[i])

      exec("pixel%d_flux = clipped_flux" % index)
      exec("pixel%d_time = clipped_time" % index)
 
      # export smoothed and clipped data as .dat file
      # exportblend = np.array([clipped_time, clipped_flux])
      # exportblend = np.transpose(exportblend)
      # exec("np.savetxt('kicskeleton_pixel%d_lc.dat', exportblend, delimiter=' ', header='Smoothed and clipped light curve for KICskeleton TPF')" % index)
      
      # fourier transform
      frequencies, power_spectrum = fft(np.asarray(clipped_time), np.asarray(clipped_flux), ofac, hifac)
      if frequencies.shape != power_spectrum.shape:
         power_spectrum = np.concatenate((power_spectrum, [0]))
      power_spectrum = np.sqrt(power_spectrum)
      power_spectrum = 1e6 * ((power_spectrum / (np.median(power_spectrum))) - 1.) # to ppm

      exec("pixel%d_freq = frequencies" % index)
      exec("pixel%d_ps = power_spectrum" % index)


### PLOTTING ###

# light curves
fig = plt.figure(1)
gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
plt.title('skeleton')
plt.xlabel('Time (d)')
plt.ylabel('Fractional Intensity')

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
   else:
      exec("flux = pixel%d_flux" % index)
      exec("time = pixel%d_time" % index)
      ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
      if img == np.amax(table2):
         plt.plot(time, flux, 'r-')
         plt.ylim(min(flux), max(flux)) #ymin=0)
         plt.xlim(min(time), max(time))
      else:
         plt.plot(time, flux, 'k-')
         plt.ylim(min(flux), max(flux)) #ymin=0)
         plt.xlim(min(time), max(time))

plt.savefig('kicskeleton_pixelslc.png')

# power spectra
fig = plt.figure(2)
gs = gridspec.GridSpec(y, x, wspace=0, hspace=0)
plt.title('skeleton')
plt.xlabel('Frequency ($\mu$Hz)')
plt.ylabel('Amplitude (ppm)')

for (j, k), img in np.ndenumerate(table2):
   index = (j + 1) * k + (x - k) * j
   if img == 0:
      ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
   else:
      exec("freq = pixel%d_freq" % index)
      exec("ps = pixel%d_ps" % index)
      ax = fig.add_subplot(gs[y - j - 1, k])
      ax.set_xticklabels('')
      ax.set_yticklabels('')
      if img == np.amax(table2):
         plt.plot(freq, ps, 'r-')
         plt.ylim(0, max(ps)) #ymin=0)
         plt.xlim(0, max(freq))
      else:
         plt.plot(freq, ps, 'k-')
         plt.ylim(0, max(ps)) #ymin=0)
         plt.xlim(0, max(freq))

plt.savefig('kicskeleton_pixels.png')
