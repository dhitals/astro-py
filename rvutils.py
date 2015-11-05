""" Utilities needed for measuring RV of MS and WD stars.

    Written by: Saurav Dhital (October 2015)
"""
import numpy as np

def readSpec(fname):
    """ read a 1D IRAF-styped spectrum
        should include fluxing at some point
        """
    from astropy.io import fits
    
    hdulist = fits.open(fname)
    hdu = hdulist[0]
    hdr = hdu.header
    
    flux = hdu.data
    
    # construct the wavelength array
    lam0 = hdr['CRVAL1']
    dlam = hdr['CD1_1'] # or CD1_1
    
    wave = lam0 + dlam *(np.arange(hdr['NAXIS1']) + 1.0 - hdr['CRPIX1'])

    return (wave, flux, hdr)


def loglinear(wave, flux):
    c = 2.99792458e5 # speed of light in km/s
    
    R = wave[0] / (wave[1] - wave[0]) ## R = lambda / dlambda -- set to the lowest value in the cc_range
    dv = c / R                     ## velocity resolution per pixel in km/s         

    ## create the new wavelength array -- log-linear
    ## Why log-linear? dv/c = dlambda/lambda = (d/dl) log_lambda
    step = (dv/c) * np.log10(np.exp(1))
    lam0, lam1 = np.log10(wave[0]), np.log10(wave[-1])
    nsteps = np.ceil( (lam1 - lam0) / step )  

    logWave = np.logspace(lam0, lam1, nsteps)
    
    logFlux = np.interp(logWave, wave, flux)
    
    return (logWave, logFlux)

# calculat the heliocentric velocity
def wave_helio(wave, hdr):
    from astropysics.coords import AngularCoordinate
    from astropy.time import Time
    from PyAstronomy import pyasl

    # coords for object -- current epoch
    ra  = AngularCoordinate(hdr['RA'],sghms=True).degrees # read RA/DEC in degrees
    dec = AngularCoordinate(hdr['DEC'],sghms=False).degrees

    # get observatory coords and JD of observation
    # use different ones for each observatory b/c their headers are different
    # currently JD is start of exposure instead of middle --- worth fixing
    if hdr['OBSERVAT'] == 'KPNO':
        lon, lat, alt, tz = 111.600, +31.9633, 2120, 7
        jd = np.double(hdr['JD'])
    elif hdr['OBSERVAT'] == 'SOAR':
        lon, lat, alt, tz = 70.7336, -30.2379, 2713., 4
        t = Time(hdr['DATE-OBS'], format='isot', scale='utc')
        jd = t.jd
    elif hdr['OBSERVAT'] == 'MMTO':
        lon, lat, alt, tz = 110.885, +31.6883, 2600, 7
        date = hdr['DATE-OBS']+'T'+hdr['UT']    # Time of observation
        t = Time(date, format='isot', scale='utc')
        jd = t.jd
    
    # Calculate barycentric correction (debug=True shows various intermediate results)
    helio_corr, hjd = pyasl.helcorr(lon, lat, alt, ra, dec, jd, debug=False)

    #gamma = sqrt( (1. + helio_corr/c) / (1. - helio_corr/c) )
    #wave *= gamma
    # print "Barycentric correction [km/s]: ", helio_corr
    #print "Heliocentric Julian day: ", hjd

    return helio_corr

def deg2rad(degrees):
    return degrees*pi/180.

def rad2deg(radians):
    return radians*180./pi

def readtem(spT):
    ## Read the template spectrum given a spectral type 

    from astropy.table import Table
    from astropy.io import fits

    ML_PATH  = '/Users/saurav/tara/local/datalib/spectral_templates/bochanski/'
    FGK_PATH = '/Users/saurav/Dropbox/wkd/grs/templates/L/'

    if spT[0].lower() == 'l':
        temName = ML_PATH + 'L0.all.na.k.fits'
        wave, flux, hdr = readSpec(temName)
    elif spT[0].lower() == 'm':
        temName = ML_PATH + spT[0:2].lower() + '.all.na.k.fits'
        wave, flux, hdr = readSpec(temName)
    else:
        temName = FGK_PATH + spT.upper() + '.fits'
        tem = np.array(fits.getdata(temName)) # elodie templates were written sans headers
        wave, flux = tem[0,:], tem[1,:]
    
    return (wave, flux)