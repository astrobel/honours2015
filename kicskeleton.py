import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits as pyfits
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel
import smoothing

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{helvet}',
       r'\usepackage[EULERGREEK]{sansmath}',
       r'\sansmath'
]


print ' '
smoothtype = raw_input('Boxcar or Gaussian smoothing? (b/g) ')
while smoothtype != "b" and smoothtype != "g":
   smoothtype = raw_input('Please type either b or g: ')
kern = raw_input('Smoothing kernel: (integer) ')
kern = np.float64(kern)

### ZEROTH QUARTER ### quarters can be removed as needed

hdulist = pyfits.open('kplr00skeleton-2009131105131_llc.fits')

table_0 = hdulist[1].data
sap_flux_0 = table_0['SAP_FLUX']
#sap_flux_err_0 = table_0['SAP_FLUX_ERR']
time_0 = table_0['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_0 = np.array([time_0, sap_flux_0])
blend_0 = np.transpose(blend_0)
blend2_0 = np.ma.compress_rows(np.ma.fix_invalid(blend_0))

time_0 = blend2_0[:,0]
sap_flux_0 = blend2_0[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_0, smth_flux_0 = smoothing.boxsmooth(time_0, sap_flux_0, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_0, smth_flux_0 = smoothing.gausssmooth(time_0, sap_flux_0, kern)


### FIRST QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2009166043257_llc.fits')

table_1 = hdulist[1].data
sap_flux_1 = table_1['SAP_FLUX']
#sap_flux_err_1 = table_1['SAP_FLUX_ERR']
time_1 = table_1['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_1 = np.array([time_1, sap_flux_1])
blend_1 = np.transpose(blend_1)
blend2_1 = np.ma.compress_rows(np.ma.fix_invalid(blend_1))

time_1 = blend2_1[:,0]
sap_flux_1 = blend2_1[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_1, smth_flux_1 = smoothing.boxsmooth(time_1, sap_flux_1, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_1, smth_flux_1 = smoothing.gausssmooth(time_1, sap_flux_1, kern)


### SECOND QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2009259160929_llc.fits')

table_2 = hdulist[1].data
sap_flux_2 = table_2['SAP_FLUX']
#sap_flux_err_2 = table_2['SAP_FLUX_ERR']
time_2 = table_2['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_2 = np.array([time_2, sap_flux_2])
blend_2 = np.transpose(blend_2)
blend2_2 = np.ma.compress_rows(np.ma.fix_invalid(blend_2))

time_2 = blend2_2[:,0]
sap_flux_2 = blend2_2[:,1]
      
if smoothtype == "b": # boxcar smoothing
   sap_flux2_2, smth_flux_2 = smoothing.boxsmooth(time_2, sap_flux_2, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_2, smth_flux_2 = smoothing.gausssmooth(time_2, sap_flux_2, kern)


### THIRD QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2009350155506_llc.fits')

table_3 = hdulist[1].data
sap_flux_3 = table_3['SAP_FLUX']
#sap_flux_err_3 = table_3['SAP_FLUX_ERR']
time_3 = table_3['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_3 = np.array([time_3, sap_flux_3])
blend_3 = np.transpose(blend_3)
blend2_3 = np.ma.compress_rows(np.ma.fix_invalid(blend_3))

time_3 = blend2_3[:,0]
sap_flux_3 = blend2_3[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_3, smth_flux_3 = smoothing.boxsmooth(time_3, sap_flux_3, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_3, smth_flux_3 = smoothing.gausssmooth(time_3, sap_flux_3, kern)


### FOURTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2010078095331_llc.fits')

table_4 = hdulist[1].data
sap_flux_4 = table_4['SAP_FLUX']
#sap_flux_err_4 = table_4['SAP_FLUX_ERR']
time_4 = table_4['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_4 = np.array([time_4, sap_flux_4])
blend_4 = np.transpose(blend_4)
blend2_4 = np.ma.compress_rows(np.ma.fix_invalid(blend_4))

time_4 = blend2_4[:,0]
sap_flux_4 = blend2_4[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_4, smth_flux_4 = smoothing.boxsmooth(time_4, sap_flux_4, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_4, smth_flux_4 = smoothing.gausssmooth(time_4, sap_flux_4, kern)


### FIFTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2010174085026_llc.fits')

table_5 = hdulist[1].data
sap_flux_5 = table_5['SAP_FLUX']
#sap_flux_err_5 = table_5['SAP_FLUX_ERR']
time_5 = table_5['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_5 = np.array([time_5, sap_flux_5])
blend_5 = np.transpose(blend_5)
blend2_5 = np.ma.compress_rows(np.ma.fix_invalid(blend_5))

time_5 = blend2_5[:,0]
sap_flux_5 = blend2_5[:,1]
 
if smoothtype == "b": # boxcar smoothing
   sap_flux2_5, smth_flux_5 = smoothing.boxsmooth(time_5, sap_flux_5, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_5, smth_flux_5 = smoothing.gausssmooth(time_5, sap_flux_5, kern)


### SIXTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2010265121752_llc.fits')

table_6 = hdulist[1].data
sap_flux_6 = table_6['SAP_FLUX']
#sap_flux_err_6 = table_6['SAP_FLUX_ERR']
time_6 = table_6['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_6 = np.array([time_6, sap_flux_6])
blend_6 = np.transpose(blend_6)
blend2_6 = np.ma.compress_rows(np.ma.fix_invalid(blend_6))

time_6 = blend2_6[:,0]
sap_flux_6 = blend2_6[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_6, smth_flux_6 = smoothing.boxsmooth(time_6, sap_flux_6, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_6, smth_flux_6 = smoothing.gausssmooth(time_6, sap_flux_6, kern)


### SEVENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2010355172524_llc.fits')

table_7 = hdulist[1].data
sap_flux_7 = table_7['SAP_FLUX']
#sap_flux_err_7 = table_7['SAP_FLUX_ERR']
time_7 = table_7['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_7 = np.array([time_7, sap_flux_7])
blend_7 = np.transpose(blend_7)
blend2_7 = np.ma.compress_rows(np.ma.fix_invalid(blend_7))

time_7 = blend2_7[:,0]
sap_flux_7 = blend2_7[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_7, smth_flux_7 = smoothing.boxsmooth(time_7, sap_flux_7, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_7, smth_flux_7 = smoothing.gausssmooth(time_7, sap_flux_7, kern)


### EIGHTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2011073133259_llc.fits')

table_8 = hdulist[1].data
sap_flux_8 = table_8['SAP_FLUX']
#sap_flux_err_8 = table_8['SAP_FLUX_ERR']
time_8 = table_8['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_8 = np.array([time_8, sap_flux_8])
blend_8 = np.transpose(blend_8)
blend2_8 = np.ma.compress_rows(np.ma.fix_invalid(blend_8))

time_8 = blend2_8[:,0]
sap_flux_8 = blend2_8[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_8, smth_flux_8 = smoothing.boxsmooth(time_8, sap_flux_8, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_8, smth_flux_8 = smoothing.gausssmooth(time_8, sap_flux_8, kern)


### NINTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2011177032512_llc.fits')

table_9 = hdulist[1].data
sap_flux_9 = table_9['SAP_FLUX']
#sap_flux_err_9 = table_9['SAP_FLUX_ERR']
time_9 = table_9['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_9 = np.array([time_9, sap_flux_9])
blend_9 = np.transpose(blend_9)
blend2_9 = np.ma.compress_rows(np.ma.fix_invalid(blend_9))

time_9 = blend2_9[:,0]
sap_flux_9 = blend2_9[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_9, smth_flux_9 = smoothing.boxsmooth(time_9, sap_flux_9, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_9, smth_flux_9 = smoothing.gausssmooth(time_9, sap_flux_9, kern)


### TENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2011271113734_llc.fits')

table_10 = hdulist[1].data
sap_flux_10 = table_10['SAP_FLUX']
#sap_flux_err_10 = table_10['SAP_FLUX_ERR']
time_10 = table_10['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_10 = np.array([time_10, sap_flux_10])
blend_10 = np.transpose(blend_10)
blend2_10 = np.ma.compress_rows(np.ma.fix_invalid(blend_10))

time_10 = blend2_10[:,0]
sap_flux_10 = blend2_10[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_10, smth_flux_10 = smoothing.boxsmooth(time_10, sap_flux_10, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_10, smth_flux_10 = smoothing.gausssmooth(time_10, sap_flux_10, kern)


### ELEVENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2012004120508_llc.fits')

table_11 = hdulist[1].data
sap_flux_11 = table_11['SAP_FLUX']
#sap_flux_err_11 = table_11['SAP_FLUX_ERR']
time_11 = table_11['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_11 = np.array([time_11, sap_flux_11])
blend_11 = np.transpose(blend_11)
blend2_11 = np.ma.compress_rows(np.ma.fix_invalid(blend_11))

time_11 = blend2_11[:,0]
sap_flux_11 = blend2_11[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_11, smth_flux_11 = smoothing.boxsmooth(time_11, sap_flux_11, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_11, smth_flux_11 = smoothing.gausssmooth(time_11, sap_flux_11, kern)


### TWELFTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2012088054726_llc.fits')

table_12 = hdulist[1].data
sap_flux_12 = table_12['SAP_FLUX']
#sap_flux_err_12 = table_12['SAP_FLUX_ERR']
time_12 = table_12['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_12 = np.array([time_12, sap_flux_12])
blend_12 = np.transpose(blend_12)
blend2_12 = np.ma.compress_rows(np.ma.fix_invalid(blend_12))

time_12 = blend2_12[:,0]
sap_flux_12 = blend2_12[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_12, smth_flux_12 = smoothing.boxsmooth(time_12, sap_flux_12, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_12, smth_flux_12 = smoothing.gausssmooth(time_12, sap_flux_12, kern)


### THIRTEENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2012179063303_llc.fits')

table_13 = hdulist[1].data
sap_flux_13 = table_13['SAP_FLUX']
#sap_flux_err_13 = table_13['SAP_FLUX_ERR']
time_13 = table_13['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_13 = np.array([time_13, sap_flux_13])
blend_13 = np.transpose(blend_13)
blend2_13 = np.ma.compress_rows(np.ma.fix_invalid(blend_13))

time_13 = blend2_13[:,0]
sap_flux_13 = blend2_13[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_13, smth_flux_13 = smoothing.boxsmooth(time_13, sap_flux_13, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_13, smth_flux_13 = smoothing.gausssmooth(time_13, sap_flux_13, kern)


### FOURTEENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2012277125453_llc.fits')

table_14 = hdulist[1].data
sap_flux_14 = table_14['SAP_FLUX']
#sap_flux_err_14 = table_14['SAP_FLUX_ERR']
time_14 = table_14['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_14 = np.array([time_14, sap_flux_14])
blend_14 = np.transpose(blend_14)
blend2_14 = np.ma.compress_rows(np.ma.fix_invalid(blend_14))

time_14 = blend2_14[:,0]
sap_flux_14 = blend2_14[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_14, smth_flux_14 = smoothing.boxsmooth(time_14, sap_flux_14, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_14, smth_flux_14 = smoothing.gausssmooth(time_14, sap_flux_14, kern)


### FIFTEENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2013011073258_llc.fits')

table_15 = hdulist[1].data
sap_flux_15 = table_15['SAP_FLUX']
#sap_flux_err_15 = table_15['SAP_FLUX_ERR']
time_15 = table_15['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_15 = np.array([time_15, sap_flux_15])
blend_15 = np.transpose(blend_15)
blend2_15 = np.ma.compress_rows(np.ma.fix_invalid(blend_15))

time_15 = blend2_15[:,0]
sap_flux_15 = blend2_15[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_15, smth_flux_15 = smoothing.boxsmooth(time_15, sap_flux_15, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_15, smth_flux_15 = smoothing.gausssmooth(time_15, sap_flux_15, kern)


### SIXTEENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2013098041711_llc.fits')

table_16 = hdulist[1].data
sap_flux_16 = table_16['SAP_FLUX']
#sap_flux_err_16 = table_16['SAP_FLUX_ERR']
time_16 = table_16['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_16 = np.array([time_16, sap_flux_16])
blend_16 = np.transpose(blend_16)
blend2_16 = np.ma.compress_rows(np.ma.fix_invalid(blend_16))

time_16 = blend2_16[:,0]
sap_flux_16 = blend2_16[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_16, smth_flux_16 = smoothing.boxsmooth(time_16, sap_flux_16, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_16, smth_flux_16 = smoothing.gausssmooth(time_16, sap_flux_16, kern)


### SEVENTEENTH QUARTER ###

hdulist = pyfits.open('kplr00skeleton-2013131215648_llc.fits')

table_17 = hdulist[1].data
sap_flux_17 = table_17['SAP_FLUX']
#sap_flux_err_17 = table_17['SAP_FLUX_ERR']
time_17 = table_17['TIME']

hdulist.close()

# creating a blend array to remove NaNs
blend_17 = np.array([time_17, sap_flux_17])
blend_17 = np.transpose(blend_17)
blend2_17 = np.ma.compress_rows(np.ma.fix_invalid(blend_17))

time_17 = blend2_17[:,0]
sap_flux_17 = blend2_17[:,1]

if smoothtype == "b": # boxcar smoothing
   sap_flux2_17, smth_flux_17 = smoothing.boxsmooth(time_17, sap_flux_17, kern)
elif smoothtype == "g": # gaussian smoothing
   sap_flux2_17, smth_flux_17 = smoothing.gausssmooth(time_17, sap_flux_17, kern)


### PLOTTING ###

# unsmoothed data
plt.figure(1)
plt.plot(   
   time_0, sap_flux_0, 'ro', time_0, smth_flux_0, 'c-',
   time_1, sap_flux_1, 'yo', time_1, smth_flux_1, 'c-',
   time_2, sap_flux_2, 'go', time_2, smth_flux_2, 'c-',
   time_3, sap_flux_3, 'bo', time_3, smth_flux_3, 'c-',
   time_4, sap_flux_4, 'mo', time_4, smth_flux_4, 'c-',
   time_5, sap_flux_5, 'ro', time_5, smth_flux_5, 'c-',
   time_6, sap_flux_6, 'yo', time_6, smth_flux_6, 'c-',
   time_7, sap_flux_7, 'go', time_7, smth_flux_7, 'c-',
   time_8, sap_flux_8, 'bo', time_8, smth_flux_8, 'c-',
   time_9, sap_flux_9, 'mo', time_9, smth_flux_9, 'c-',
   time_10, sap_flux_10, 'ro', time_10, smth_flux_10, 'c-',
   time_11, sap_flux_11, 'yo', time_11, smth_flux_11, 'c-',
   time_12, sap_flux_12, 'go', time_12, smth_flux_12,'c-',
   time_13, sap_flux_13, 'bo', time_13, smth_flux_13, 'c-',
   time_14, sap_flux_14, 'mo', time_14, smth_flux_14, 'c-',
   time_15, sap_flux_15, 'ro', time_15, smth_flux_15, 'c-',
   time_16, sap_flux_16, 'yo', time_16, smth_flux_16, 'c-',
   time_17, sap_flux_17, 'go', time_17, smth_flux_17, 'c-', markersize=3
   )
plt.xlabel('Time (d)')
plt.ylabel('Flux (e$^{-}$/sec)')
if smoothtype == "b":
   plt.title('skeleton: Raw with Boxcar Fit')
elif smoothtype == "g":
   plt.title('skeleton: Raw with Gaussian Fit')
plt.savefig('kicskeleton_raw.png')

# concatenation for smoothed data
time = np.concatenate((time_0, time_1, time_2, time_3, time_4, time_5, 
   time_6, time_7, time_8, time_9, time_10, time_11, time_12, time_13,
   time_14, time_15, time_16, time_17))
sap_flux = np.concatenate((sap_flux2_0, sap_flux2_1, sap_flux2_2,
   sap_flux2_3, sap_flux2_4, sap_flux2_5, sap_flux2_6, sap_flux2_7,
   sap_flux2_8, sap_flux2_9, sap_flux2_10, sap_flux2_11, sap_flux2_12,
   sap_flux2_13, sap_flux2_14, sap_flux2_15, sap_flux2_16, sap_flux2_17))

# scan for outliers
print ' '
inp = raw_input('Clipping level: ')
print ' '
clip = np.float64(inp) * np.std(sap_flux)
meanflux = np.mean(sap_flux)

upperbound = meanflux + clip
lowerbound = meanflux - clip

colours = np.zeros(sap_flux.size)

for i, flux in enumerate(sap_flux):
   if flux < upperbound and flux > lowerbound:
      colours[i] = 1

clipped_flux = []
clipped_time = []

# smoothed data
plt.figure(2)
plt.plot(time, sap_flux, 'bo', markersize=3)
plt.xlabel('Time (d)')
plt.ylabel('Fractional Intensity')
plt.title('skeleton')
for i, colour in enumerate(colours):
   if colour == 1:
      clipped_flux.append(sap_flux[i])
      clipped_time.append(time[i])
plt.plot(clipped_time, clipped_flux, 'ro', markersize=3)
plt.savefig('kicskeleton_smooth.png')

# export smoothed and clipped data as .dat file
exportblend = np.array([clipped_time, clipped_flux])
exportblend = np.transpose(exportblend)
np.savetxt('kicskeleton_lc.dat', exportblend, delimiter=' ',
   header='Smoothed and clipped light curve for KICskeleton')

# plt.show()
