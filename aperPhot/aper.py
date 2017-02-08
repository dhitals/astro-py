""" Perform aperture photometry on a set of images given x,y coords

    Make sure the input, telescope, and photomotery params are correct.
    They are all specified at the beginning of this file.

    Written by: S Dhital (June 2015)
"""

import numpy as np
import astropy.table
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.table import hstack
from tabulate import tabulate
import glob, os

## camera parameters -- for SARA ARC cameras
gain = 2.0 # e-/ADU 
rdnoise = 6.0 # ADU

# aperture photometry parameters
aper_size = 8.
annulus, dannulus = 15., 20.

zmag = 25.  # zero-point of magnitude scale 0 -- arbitrary, obviously

# all of the input filenames -- make sure that they are correctly specified
baseFname = 'PG1613+426'
path = 'images/'

coordFile = 'ref.' + baseFname + '.coords' # list of coordinates (x, y)
photFile = path + baseFname + '.phot' # output filename

# list of images to do photometry on 
imFiles = glob.glob(path + 'Sh.' + baseFname + '*.fits')
nImages = len(imFiles)
print('No. of images in this dir: %s' %(nImages))

# read in coordinates of the stars -- assumes (N,2)
stars = np.genfromtxt(coordFile, names=True, dtype='f8')
nStars = stars.size
print('No. of stars in ref image: %i\n' %(nStars))
coords = (stars['x'], stars['y'])

# read in the FWHM of the images -- output from kfwhm
fwhms = np.genfromtxt('kfwhm.dat', names=True, dtype='S30,f8,f8,f8,i8')
print tabulate(fwhms,fwhms.dtype.names)


""" Aperture photometry happens here.
    follow the algorithm outlined in 
    http://photutils.readthedocs.org/en/latest/photutils/aperture.html#global-background-subtraction
    error in photometry is assumed as Poissonian -- sqrt(N)
    but smooth background error not included for now
    """
for eachFile in imFiles:
    print('Prcoessing %s ...' %(eachFile))
    
    try: 
        hdu = pyfits.open(eachFile)
        im = hdu[0].data
    except: print('File could not be opened: %s' %(eachFile))

    # create aperture and annulus objects
    apertures = CircularAperture(coords, r=aper_size)
    annulus_apertures = CircularAnnulus(coords, r_in=annulus, r_out=dannulus)
    
    npix_src, npix_bkg = apertures.area(),  annulus_apertures.area()
    
    # calculate the object and annulus flux
    data_error = np.sqrt(im)
    rawflux_table = aperture_photometry(im, apertures, error=data_error)
    bkgflux_table = aperture_photometry(im, annulus_apertures, error=data_error)
    # stack the two tables into one
    phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])

    # calculate bkg mean & normalize by area
    bkg_mean = phot_table['aperture_sum_bkg'] / npix_bkg

    # final photometric counts -- in ADUs
    Nsrc = phot_table['aperture_sum_raw'] - bkg_mean * npix_src
    # change counts to flux
    flux = gain * Nsrc / hdu[0].header['EXPTIME']

    # calculate the photometric error (see CCD Eq in Howell's Hanbook of CCD Astronomy):
    # ignore noise due to dark current and error due to A/D converter
    # add corrections for co-add (N_i) and if median is used (k)
    #        sigma_src^2 = (S_i - <B>) / N_i + n_A * (1 + k * n_A/n_B) * (N_sky + N_rdnoise^2)
    #
    # S_i => counts per pixel
    # <B> => mean bkg count
    # N_i => # of co-adds for the pixels. Assume N_i = 1
    # n_A => No. of pixels in source aperture
    # n_B => No. of pixels in bkg aperture
    # k = 1 if <B> is mean. k = pi/2 if <B> is median
    
    fluxerr = np.sqrt(Nsrc + 
                      npix_src * (1 + npix_src/npix_bkg) * (bkg_mean + np.square(rdnoise)) )
    
    # calculate mag, magerr from counts
    phot_table['mag'] = zmag - 2.5 * np.log10(flux)
    phot_table['magerr'] = 1.0857 * fluxerr / Nsrc  ## approximation from Howell 1993 IAUS

    #print phot_table['magerr'][0],fluxerr[0], Nsrc[0], npix_src, npix_bkg, phot_table['aperture_sum_bkg'][0],np.square(rdnoise)
    
    # Now, lets add the following cols to the table for the phot.out file
    # Required for output:
    #                   imageName, starNum, JD, mag, magErr
    # Optional:
    #                   filter, dateObs, exptime, AM, fwhm 

    #print phot_table.columns
    fname = eachFile.split('/')[-1]
    photFile = fname.split('.fits')[0]+'.phot'
    
    # get this from coordfile
    starNum = astropy.table.Column(name='starNum', data=stars['starNum'])

    imageName = astropy.table.Column(name='imageName', data=np.chararray(nStars) + eachFile.split('/')[-1])
    # get FWHM    
    ind = np.where(fwhms['frame'] == path+fname.split('.')[-2]+'.fits')  ## HARD-CODED
    fwhm_star = 0.0 if len(ind[0]) == 0 else fwhms['FWHM'][ind]
    fwhm = astropy.table.Column(name='fwhm', data=np.zeros(nStars) + fwhm_star)

    # get these from the header
    jdate   = astropy.table.Column(name='jdate',   
                                   data=np.ones(nStars) +  float(hdu[0].header['J_DATE'].split()[0]))
    obsdate = astropy.table.Column(name='obsdate', data=np.chararray(nStars) + hdu[0].header['DATE-OBS'])
    airmass = astropy.table.Column(name='airmass', data=np.ones(nStars) + hdu[0].header['AIRMASS'])

    # add using add_columns s.t. the fields are in proper order
    phot_table.add_columns([imageName,jdate, airmass, fwhm, starNum], [0,0,0,0,0])
    
    # rename these cols
    phot_table.rename_column('xcenter_raw', 'xcenter')
    phot_table.rename_column('ycenter_raw', 'ycenter')

    # don't print these cols
    phot_table.remove_columns(('xcenter_bkg', 'ycenter_bkg', 'aperture_sum_bkg','aperture_sum_err_bkg',
                               'aperture_sum_raw','aperture_sum_err_raw'))

    #print phot_table.columns
    #print tabulate(phot_table, headers=phot_table.colnames)
    np.savetxt(path+photFile, phot_table, #header=phot_table.colnames,
               fmt='%s %.5f %5.2f %5.2f %03d %7.1f %7.1f %8.2f %8.2f')
    print '         ... Done'
    
os.system("echo 'imageName             jdate       airmass fwhm starNum xcenter ycenter  mag    magerr' > phot.out") 
os.system("cat images/Sh*.phot >> phot.dat")  ## PATH hard-coded -- FIX
print "Phot file written to phot.out"

