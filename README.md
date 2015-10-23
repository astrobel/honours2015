# honours2015
This is a collection of the code I have used for my Physics IV 2015 research project. It is written in Python 2.7, using packages AstroPy, NumPy, SciPy and matplotlib.

Filenames and files contain the word "skeleton" (indicating that they are the bones of workable code) which is to be replaced with the KIC number of the star being analysed. Future improvements to this code will involve automation so that this does not have to be changed manually every time. Note also that input files must be added manually.

Description of code included:

  My code:

    - smoothing.py: includes routines for polynomial, boxcar, and Gaussian smoothing methods

    - kicskeleton.py: SAP light curve cleaning and smoothing

    - kicskeleton_ft.py: handles Fourier transform of cleaned light curve and produces amplitude spectrum

    - kicskeleton_pc.py: creates a zoomed-in amplitude spectrum, designed to be generated with a higher oversampling factor than in the previous code, and uses the highest peak to fold a phased light curve

    - kicskeletonpixel_even.py: analyses TPF light curves (cleaning, smoothing, Fourier analysis) for stars located on even channels of the Kepler telescope


    - kicskeletonpixel_odd.py: analyses TPF light curves for stars located on odd channels of the Kepler telescope

    - kicskeletonimg_evenpy: compares TPF and UKIRT images for stars located on even channels of the Kepler telescope


    - kicskeletonimg_odd.py: compares TPF and UKIRT images for stars located on odd channels of the Kepler telescope

    - kicskeletonloc.py: generates two plots, one of the UKIRT image of the target and one of this image with KIC and UKIRT catalogue cross-matching overlayed
    
  Code from elsewhere (links in files):

    - lomb.py: Lomb-Scargle Fourier transform routine

    - translate.py: routine to map one range of values to another

  Extras:

    - quarterstitchingpatch.txt: added to pixel analysis code in the case that multiple quarters of data were used

    - refpix.txt: added to TPF image analysis to locate the precise reference coordinates used for a target star
