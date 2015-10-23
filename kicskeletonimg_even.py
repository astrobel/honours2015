import numpy as np
from astropy.io import fits as pyfits
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{helvet}',
       r'\usepackage[EULERGREEK]{sansmath}',
       r'\sansmath'
]

### IMAGE 1: ONE KEPLER PIXEL IMAGE, Q_ ###

hdulist = pyfits.open('*', mode='readonly')

table = hdulist[1].data
flux = table['FLUX']
#flux_err = table['FLUX_ERR']
time = table['TIME']
hd1 = hdulist[1].header

w = wcs.WCS(hd1, keysel=['binary'])

table2 = hdulist[2].data
hd2 = hdulist[2].header
x = hd2['NAXIS1'] 
y = hd2['NAXIS2']

hdulist.close()

imgflux = np.flipud(flux[0])
imgflux = np.fliplr(imgflux)


### IMAGE 2: UKIRT IMAGE ###

hdulist = pyfits.open('*')

flux2 = hdulist[1].data


### PLOTTING ###

fig, (kepler, ukirt) = plt.subplots(1, 2) 

fig.suptitle('skeleton', fontsize=20)

left = kepler.imshow(imgflux, cmap='pink')
kepler.set_title('Kepler Aperture')
left.set_interpolation('nearest')
kepler.set_xlim(-0.5, x-0.5)
kepler.set_ylim(y-0.5, -0.5)

left.axes.get_xaxis().set_ticklabels([])
left.axes.get_yaxis().set_ticklabels([])
left.axes.get_xaxis().set_ticks([])
left.axes.get_yaxis().set_ticks([])

crval = w.wcs.crval
north = crval + np.array([0, 6/3600.])
east = crval + np.array([ 6/3600., 0])

ncoords = np.vstack([crval, north])
ecoords = np.vstack([crval, east])
npixels = w.wcs_world2pix(ncoords , 0)
epixels = w.wcs_world2pix(ecoords , 0)
npixels[1, 1] = npixels[0, 1] - (npixels[1, 1] - npixels[0, 1]) # flip ud
epixels[1, 1] = epixels[0, 1] - (epixels[1, 1] - epixels[0, 1]) #+ 0.2
# epixels[1, 0] += 0.2
npixels[1, 0] = npixels[0, 0] - (npixels[1, 0] - npixels[0, 0]) # flip lr
epixels[1, 0] = epixels[0, 0] - (epixels[1, 0] - epixels[0, 0])
kepler.plot(npixels[:,0], npixels[:,1], color='#00ff8c')
kepler.plot(epixels[:,0], epixels[:,1], color='#00ff8c')

kepler.text(npixels[1, 0] - 0.3, npixels[1, 1] + 0.4, 'N', color='#00ff8c')
kepler.text(epixels[1, 0], epixels[1, 1] + 0.4, 'E', color='#00ff8c')

right = ukirt.imshow(flux2, cmap='pink', norm=LogNorm())
ukirt.set_title('UKIRT Image')
right.set_interpolation('bilinear')
ukirt.set_xlim(300, 0)
ukirt.set_ylim(0, 300)

right.axes.get_xaxis().set_ticklabels([])
right.axes.get_yaxis().set_ticklabels([])
right.axes.get_xaxis().set_ticks([])
right.axes.get_yaxis().set_ticks([])

ukirt.plot([25, 25], [25, 55], '-', color='#00ff8c')
ukirt.plot([25, 55], [25, 25], '-', color='#00ff8c')

ukirt.text(30, 60, 'N', color='#00ff8c')
ukirt.text(70, 20, 'E', color='#00ff8c')


fig.set_size_inches(7.5, 4.5)
plt.savefig('kicskeletonimg.png')

# plt.show()
