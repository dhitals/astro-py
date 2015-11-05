""" Use cross-correlation to measure RV of main-sequence stars.
    Use a sliding box to measure the shift and averages the ones that converge
    Some bugs still need to be worked out
    
    Written by Saurav Dhital (November 2015)
    """

from __future__ import print_function
import glob
import numpy as np
import scipy as sp
import atpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib import gridspec
get_ipython().magic(u'matplotlib inline')

from rvutils import *
import seaborn as sns

def set_axes():

    fig = plt.figure(figsize=(12, 6), dpi=200) 

    # split the figure vertically in half
    outer_grid = gridspec.GridSpec(1, 2, wspace=0.1, width_ratios=[1,1])

    # the LHS is for the spectrum
    ax0 = plt.subplot(outer_grid[0])
    ax0.set_xlabel(r'$\lambda (\AA$)')
    ax0.set_ylabel('Flux')
    
    # grid the RHS for CCF
    nx, ny = 3, 3
    inner_grid = gridspec.GridSpecFromSubplotSpec(nx, ny, subplot_spec=outer_grid[1], wspace=0.0, hspace=0.0)

    for i in range(nx*ny):
        ax = plt.Subplot(fig, inner_grid[i])

        ax.set_xticks([-10,0,10])
        ax.set_yticks([])
        ax.set_ylim([0,10])
        ax.set_xlim([-15,15])
        ax.set_xlabel('# pixels')
        
        fig.add_subplot(ax)

    all_axes = fig.get_axes()
    
    return all_axes


def measure_rv(wdnum, spT, files):
    from PyAstronomy import pyasl    
    
    ## define paramters for cross-correlation
    c = 2.99792458e5                      ## speed of light in km/s
    box = 10.                             ## size of moving "box" for correlation
    bbox = 10*box                         ## size of window to CC
    ccRange = np.arange(5850.,6750.,bbox/2.) ## start of CC region -- 5850:6750
    nCC = np.size(ccRange)
    maxRV = 200 ## maximum possible value of RV -- absolute value

    # read the list of the spectra for that WD
    #wdnum = 'WD2350-081'
    #spT = 'K5V'

    # read the template spectra appropriate for that MS
    temWave, temFlux = readtem(spT)

    for fname in files:                    
        obsdate = fname[fname.find('_')+1 : fname.find('.fits')]
        wave, flux, hdr = readSpec(fname)    # read the object spectrum
        
        helio_corr = wave_helio(wave, hdr)   # correction to heliocentric velocity

        for j in np.arange(nCC-1):
            
            ind = (wave > ccRange[j]) & (wave < (ccRange[j] + bbox)) # select wave range

            if wave[ind].size > 1:
                c = 2.99792458e5 # speed of light in km/s
                R = wave[ind][0] / (wave[ind][1] - wave[ind][0])
                dv = c / R      ## velocity resolution per pixel in km/s         

                rv, cc = pyasl.crosscorrRV(wave[ind], flux[ind], temWave, temFlux, 
                                           -1.*maxRV, maxRV, dv, mode='doppler')

                # Find the index of maximum cross-correlation function
                maxind = np.argmax(cc)

                if abs(rv[maxind]) < maxRV:
                    print('%12s %3s %10s %.2f %i--%i %.2f %.2f' %(wdnum, spT, obsdate, helio_corr, ccRange[j], (ccRange[j]+bbox), 
                            dv, rv[maxind]+helio_corr))
                    #plt.plot(rv+helio_corr, cc, 'bp-')
                    #plt.plot(rv[maxind]+helio_corr, cc[maxind], 'ro')

                #set_axes()


if __name__ == "__main__":

    wdms = atpy.Table('../grs.fits')
    
    for i, wdnum in enumerate(wdms.wdnum):
        wdnum = wdnum.replace(' ','')
        spT = wdms['spType_ms'][i][:2]+'V'
        
        files = glob.glob('../sp-wdms/' + wdnum + 'B*.fits')
        if len(files) > 0 and spT != '\NV':
            
            ## for now, lets build in some exception for files w/ problems
            ## first problem is to get *proper* spT.
            if wdnum == 'WD0625+415': files = files[1:] # ignore the first file -- it is missing the first ~200 pixels
            if wdnum == 'WD0628-020': spT = 'M0V' # spT in grs,fits was MVV
            if wdnum == 'WD1204+450': files = files[0:2]

            print( "%2i %12s %5s N = %i" %(i, wdnum, spT, len(files)) )

            measure_rv(wdnum.replace(' ',''), spT, files)



# average over multiple CC regions to get one Rv per star

wdms = np.genfromtxt('rvMS.dat', dtype='S12,S3,S8,f8,S10,f8,f8',names=True)
print wdms.dtype.names


for wd in np.unique(wdms['wdnum']):
    ind1 = np.where(wdms['wdnum'] == wd)
    
    # now, for each obsdate, average the measured RVs
    for obsdate in np.unique(wdms[ind1]['obsdate']):
        ind2 = np.where(wdms[ind1]['obsdate'] == obsdate)
        print('%12s %8s %3i %8.2f %8.2f %8.2f' %( wd, obsdate, wdms[ind1][ind2].size, 
                                               np.mean(wdms[ind1][ind2]['rv']), 
                                               np.median(wdms[ind1][ind2]['rv']), 
                                               np.std(wdms[ind1][ind2]['rv']) ))
    
    print '#'
