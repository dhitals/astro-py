""" Measure the radial velocity of WDs using the algorithm in Becker et al. (2015)
    Much more robust than XC, which is much more robust that any other algorithm 
    out there for WDs. Hard since there is ONE single Halpha line.
    
    Yields RVs w/ posteriot dist. stddev of 5-10 km/s. 
    Mostly consistent across different epochs.
    
    Written by: Saurav Dhital, Dan (June 2015)
"""
import glob
import numpy as np
import matplotlib.pyplot as pl
import atpy
import emcee
from emcee.utils import sample_ball
from astropysics.coords import AngularCoordinate
from PyAstronomy import pyasl
from astropy.io import fits

c = 2.99792458e5 # speed of light in km/s
niter = 1e3  # Number of mcmc iterations
nwalk = 256  # Number of emcee `walkers'
gwave = (6500, 6600) # wave range to run the fit over

def deg2rad(degrees):
    return degrees*pi/180.

def rad2deg(radians):
    return radians*180./pi

def wave_helio(wave, hdr):
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

def readSpec(fname):
    """ read a 1D IRAF-styped spectrum
        return heliocentric wavelength
        should include fluxing at some point
        """
    hdulist = fits.open(fname)
    hdu = hdulist[0]
    hdr = hdu.header
    
    flux = hdu.data
    
    # construct the wavelength array
    wave = hdr['CRVAL1'] + hdr['CDELT1']*(np.arange(hdr['NAXIS1']) + 1.0 - hdr['CRPIX1'])

    # does log-linear wave help for this method?
    logWave = loglinear (wave)
    logFlux = np.interp(logWave, wave, flux)

    return (logWave, logFlux, hdr)

def loglinear(lam):
    R = lam[0] / (lam[1] - lam[0]) ## R = lambda / dlambda -- set to the lowest value in the cc_range
    dv = c / R                     ## velocity resolution per pixel in km/s         

    ## Why log-linear? dv/c = dlambda/lambda = (d/dl) log_lambda
    step = (dv/c) * np.log10(np.exp(1))
    nsteps = np.ceil( (np.log10(gwave[1]) - np.log10(gwave[0])) / step )  
    
    ## create the new wavelength array -- log-linear
    logLam = np.logspace(np.log10(gwave[0]), np.log10(gwave[1]), nsteps)
    
    return logLam    

class WDModel(object):

    def __init__(self, wave, spec):
        """Takes as arguments the model spectrum as wavelength and flux arrays
        """
        self.wave = wave
        self.flux = spec
        
    def calibration(self, norm, coeffs, wave, **kwargs):
        """
        Implements a polynomial calibration model.

        :returns cal:
           a polynomial of the form
           norm * (1 + \Sum_{i=1} * coeffs[i-1] * x **i)
        """        
        # Should find a way to make this more generic,
        # using chebyshev polynomials
        x = (wave - wave.min())/(wave.max() - wave.min())
        poly = np.zeros_like(x)
        powers = np.arange( len(coeffs) ) + 1
        poly = (x[None,:] ** powers[:,None] *
                    coeffs[:,None]).sum(axis = 0)
            #switch to have spec_norm be multiplicative or additive
            #depending on whether the calibration model is
            #multiplicative in exp^poly or just poly, Should move this
            #to mean_model()?
        return (1.0 + poly) * norm
        
    def spec(self, rv, wave):
        """Redshift the model spectrum and interpolate onto the
        observed wavelength grid
        """
        a = 1 + rv/c
        flux = np.interp(wave, self.wave * a, self.flux)
        return flux

def lnprob(theta, obs=None, model=None):
    """calculate the posterior probability
    """
    rv, norm, snr = theta[0], theta[1], theta[2]
    coeffs = theta[3:]

    # get the prior probability
    lnp_prior = prior(rv, norm, snr, coeffs)
    if not np.isfinite(lnp_prior):
        return -np.inf
    # get the model
    spec = model.spec(rv, obs['wave']) * model.calibration(norm, coeffs, obs['wave'])
    # compute the residual
    delta = spec - obs['flux']
    # estimate the uncertainty
    unc = obs['flux'] / snr
    # calculate the likelihood
    lnp = -0.5 * ((delta/unc)**2).sum() - 0.5 * np.log(unc**2).sum()
    
    return lnp + lnp_prior

def prior(rv, norm, snr, coeffs):
    """Prior probabilities for the parameters.  Returns -inf if any
    parameter is outside the priors
    
    """
    bad = ((rv < -200) | (rv > 200) |
           (norm < 0)  |
           (snr < 10)  | (snr > 100)
           )
    if bad:
         return -np.inf
    else:
        return 0

def measure_rv(wdnum):
    # center of the initial parameter guesses.  The parameters are
    # [ rv, normalization, SNR, polynomial_coeffs...., ]
    p0 = np.array([100.0, 1e19, 25.0, 1e-2, 1e-2, 1e-2])
    initial = sample_ball(p0, p0 * 0.4, nwalk)
    
    ndim = len(p0)

    # read the list of the spectra for that WD
    files = glob.glob('../sp-wdms/' + wdnum + 'A*.fits')
    print( "N = %i" %(len(files)) )
    
    for fname in files:

        # Dictionary of observational data
        obs = {}
        wave, flux, hdr = readSpec(fname)
        # Restrict the wavelength range
        gw = (wave > gwave[0]) & (wave < gwave[1])
        obs['wave'] = wave[gw]
        obs['flux'] = flux[gw]

        # model spectrum
        modelFname = '../sp-wdModel/'+wdnum+'A.model.txt'
        md = np.loadtxt(modelFname, skiprows=1)

        # instantiate the model object
        model = WDModel(md[:,0], md[:,1])
        
        # pass the model object and obs dictionary to the probability function
        lnp_kwargs = {'model':model, 'obs':obs}
        sampler = emcee.EnsembleSampler(nwalk, ndim, lnprob, kwargs=lnp_kwargs)

        # actually run the sampler
        # ERROR: NaN value of lnprob for parameters: 
        pos, prob, state = sampler.run_mcmc(initial, niter)

        pl.hist(sampler.flatchain[:,0], range=(-100, 100), bins=100)
    
        # actually render the plot
        fig, axes = pl.subplots(2,3)
        for j, p in enumerate(p0):
            ax = axes.flat[j]
            for i in range(nwalk):
                ax.plot(sampler.chain[i,:,j])
        
        # correct to heliocentric frame
        helio_corr = wave_helio(wave, hdr)

        rvArray = sampler.chain[:,500:,0]
        snrArray = sampler.chain[:,500:,2]

        t = np.array( [np.mean(rvArray)  + helio_corr,
                       np.median(rvArray)+ helio_corr,
                       np.std(rvArray),
                       np.mean(snrArray),
                       np.median(snrArray),
                       np.std(snrArray)] )

        names = ('fname','meanRV', 'medianRV', 'stdRV',
                 'helioCorr','meanSNR', 'medianSNR', 'stdSNR')
        print( "%s %8.2f %8.2f %5.2f %8.2f %10.2f %5.2f %5.2f"               %(fname, t[0], t[1], t[2], helio_corr, t[3], t[4], t[5]) )

        # save the emcee plot
        x = fname.find('WD')
        fout = fname[x:].replace('fits','png')
        pl.savefig(fout,format='png',dpi=100)

###############################################################################################

if __name__ == "__main__":

    wdms = atpy.Table('../grs.fits')
    print wdms.shape
    
    #measure_rv('WD0810+234')
    for i, wdnum in enumerate(wdms.wdnum):
        if (i > 55): #wdnum.replace(' ','') != 'WD0810+234':
            print ''
            print i, wdnum
            measure_rv(wdnum.replace(' ',''))
