#!/opt/local/bin/python
import numpy as np
import sys, os, os.path, gc, glob, datetime
from astropy.table import Table
import astropy.stats as stats
#from astropy.io import ascii
import astropy.io.fits as fits
import matplotlib.pyplot as pl
from matplotlib.colors import LogNorm

""" 
    -- Routines to do basic bias subtraction and flat fielding for a series of images
    -- Assumes that all images have been inspected for quality before running this program
    -- I recommend dokfwhm.pl, statclipped for the inspections
    -- Please manually check that reduced images are not distorted
    -- Currently bas pixel masks not incorporated
    
    -- Currently, biases and flats are not identified based on the header but by their name,
        specific to SARA images taken in 2015
"""

### paths for input & output images
inPath = '/Users/saurav/data/saraN/2015/20150611/'
outPath = '/Users/saurav/estrella/wkd/phot/L101-16/'
username = 'saurav'

obj = 'LP101-16'
biasString = 'Bias*.fits'
flatString = 'Flat*.fits'
objString = obj + '*.fits'

biasName = 'bias.fits'
flatName = 'flat.fits'

filters = ['None']

### usuaully you won't need to modify below here.


# Tweak params for source extraction if we need to
#thres = 10.0   # 5. usually  -- max allowed FWHM in pixels
#fwhm  =  3.0   # 3.0 usually

# Telescope and cameras-specific parameter
""" for SARA ARC cameras. Bill Keel is not sure if the units are photons or ADU """
#gain = 2.0 # e/ADU
#rdnoise = 6.0 # ADU

# Try to make a directory for the reduced data
#reducedPath = path+date+'/Reduced_Data/'
#try: os.mkdir(reducedPath)
#except: pass


""" print stats for each image --- name, npix, mean std, min, max """
def print_stats(im, fname):
    print('%s \t %i10 \t %0.2f \t %0.2f \t %+i \t %i'
          %(fname, im.size, np.mean(im), np.std(im), np.min(im), np.max(im) ))
    return

""" trim the overscan region -- typically specified in headers but not for SARA """
def trim_overscan(im):
    ## when specified read overscan region form header

    ### for SARA-S
    #overscan = [1000,1024]

    return im[0:1024,0:1000]


""" lets make a master bias image """
def biascombine():
    # read in biases
    files = glob.glob(inPath + biasString)
    
    biasList = []
    print '\n\nImage\t\t Npix \t\t mean \t\t stddev \t min \t max'
    
    for eachFile in files:
        try: im = np.array( fits.getdata(eachFile, ignore_missing_end=1) )
        except: print('File not found: %s' %(eachFile))
        
        im = trim_overscan(im)
        
        print_stats(im, eachFile.split('/')[-1]) # print the vitals of the image to terminal
        
        biasList.append(im)
    
    biasList = np.array(biasList)  # convert to numpy array
    
    bias = np.median(biasList, axis=0)   # median add the biases

    fits.writeto(outPath+biasName, bias, clobber=True) # write the bias
    print('Bias written to %s ' %(biasName))
    
    return bias


""" combine the flats """
def flatcombine(bias, badpix=None):
    # read in flats
    files = glob.glob(inPath + flatString)
    
    nflat = len(files)

    flatList = []
    print '\n\nImage\t\t\t Npix \t\t mean \t\t stddev \t min \t max'
    
    for eachFile in files:
        try: im = np.array( fits.getdata(eachFile, ignore_missing_end=1) )
        except: print('File not found: %s' %(eachFile))

        im = trim_overscan(im)
        
        # subtract the bias
        im -= bias

        im[im <= 0] = 1e-6 # we don't want negative values. maybe need a better fix
        
        print_stats(im, eachFile.split('/')[-1]) # print the vitals of the image to terminal

        
        flatList.append(im)
    
    flatList = np.array(flatList)

    flat = np.median(flatList, axis=0)

    # normalize the flat
    normFlat = flat / np.median(flat.flatten())

    # Write the master
    fits.writeto(outPath+flatName, normFlat, clobber=True)

    print('Flat written to %s ' %flatName)
    return normFlat


""" bias-subtract and flat-field all the images """
def ccdproc(bias, flat):

    files = glob.glob(inPath + objString) # list of all files in the dir

    for eachFile in files:
        
        try: hdr = fits.getheader(eachFile, ignore_missing_end=1)
        except:
            print('Header not readable for: %s' %(eachFile))
            continue
        
        try:
            im1 = np.array( fits.getdata(eachFile, ignore_missing_end=1) )
            im1 = trim_overscan(im1)

            im = (im1 - bias) / flat

            # update the header files with what you have done
            hdr['COMMENT'] = ( '   Bias subraction done w/ image %s' %(biasName) )
            hdr['COMMENT'] = ( '   Flat fielding   done w/ image %s' %(flatName) )
            hdr['COMMENT'] = ( '   Reduced with reduce.py on %s by %s' %(datetime.datetime.now().isoformat(), username) )
            
            # Write the master
            fits.writeto(outPath+eachFile.split('/')[-1], im, hdr, clobber=True)

        except: print('File not found: %s' %(eachFile))

if __name__ == "__main__":
    
    bias = biascombine()
    flat = flatcombine(bias)

    ccdproc(bias, flat)