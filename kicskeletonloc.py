import numpy as np
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import translate as tr

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
       r'\usepackage{helvet}',
       r'\usepackage[EULERGREEK]{sansmath}',
       r'\sansmath'
]


### UKIRT IMAGE ###

hdulist = pyfits.open('*')

flux2 = hdulist[1].data
flux2 = np.rot90(flux2, 1)


### READING IN SOME OTHER FILES ###
uk = np.loadtxt('ucoords.dat')
kep = np.loadtxt('kcoords.dat')

kepra = kep[:,0]
kepdec = kep[:,1]
kepmag = kep[:,2]
ukra = uk[:,0]
ukdec = uk[:,1]
ukmag =  uk[:,2]

for i, val in enumerate(kepmag):
   j = ukmag[i]
   if val <= 16.7:
      kepmag[i] = j - 398.04666 + 149.08127*j - 21.952130*(j**2) + 1.5968619*(j**3) - 0.057478947*(j**4) + 0.00082033223*(j**5)
   elif val > 16.7:
      kepmag[i] = j + 0.1918 + 0.08156*j


### PLOTTING ###

plt.figure(1) 

ukirt = plt.imshow(flux2, cmap='pink', norm=LogNorm())
plt.title('skeleton')
ukirt.set_interpolation('bilinear')
plt.xlim(300, 0)
plt.ylim(0, 300)

ukirt.axes.get_xaxis().set_ticklabels([])
ukirt.axes.get_yaxis().set_ticklabels([])
ukirt.axes.get_xaxis().set_ticks([])
ukirt.axes.get_yaxis().set_ticks([])

plt.plot([25, 25], [25, 55], '-', color='#00ff8c')
plt.plot([25, 55], [25, 25], '-', color='#00ff8c')

plt.text(30, 60, 'N', color='#00ff8c')
plt.text(70, 20, 'E', color='#00ff8c')

main = plt.Circle((149, 149), 15, color='#00ff8c', fill=False)
plt.gca().add_artist(main)
# opt_ = plt.Circle((x, y), 8, color='#00ff8c', fill=False)
# plt.gca().add_artist(opt_)
# plt.text(x-12, y, 'KIC', color='#00ff8c')

plt.savefig('kicskeletonloc.png')

cutoffra = 0.0107
cutoffdec = 0.0166667 / 2 # 1 arcmin in degrees to same sig fig as below... divided by 2
centra = 287.067369 #CHANGE
centdec = 39.66424

clipra1 = []
clipdec1 = []
clipcol1 = []

clipra2 = []
clipdec2 = []
clipcol2 = []

for i, ra in enumerate(ukra):
    for j, dec in enumerate(ukdec):
        if i == j and ra < centra + cutoffra and ra > centra - cutoffra and dec < centdec + cutoffdec and dec > centdec - cutoffdec:
            clipra1.append(ra)
            clipdec1.append(dec)
            clipcol1.append(ukmag[i])

for i, ra in enumerate(kepra):
    for j, dec in enumerate(kepdec):
        if i == j and ra < centra + cutoffra and ra > centra - cutoffra and dec < centdec + cutoffdec and dec > centdec - cutoffdec:
            clipra2.append(ra)
            clipdec2.append(dec)
            clipcol2.append(kepmag[i])

for i, val in enumerate(clipra1):
    clipra1[i] = tr.translate(val, centra-cutoffra, centra+cutoffra, 0, 300)
for i, val in enumerate(clipra2):
    clipra2[i] = tr.translate(val, centra-cutoffra, centra+cutoffra, 0, 300)
for i, val in enumerate(clipdec1):
    clipdec1[i] = tr.translate(val, centdec-cutoffdec, centdec+cutoffdec, 0, 300)
for i, val in enumerate(clipdec2):
    clipdec2[i] = tr.translate(val, centdec-cutoffdec, centdec+cutoffdec, 0, 300)

plt.scatter(clipra1, clipdec1, s=50, c=clipcol1, cmap='gist_rainbow', linewidths=0, alpha = 0.5)
plt.scatter(clipra2, clipdec2, marker=u'*', s=60, c=clipcol2, cmap='gist_rainbow', linewidths=0.5, edgecolor='white')#, alpha = 0.7)

cbar = plt.colorbar()
cbar.ax.invert_yaxis()

plt.savefig('kicskeletonloc1.png')

# plt.show()
