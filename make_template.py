import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import pysynphot as psyn
from astropy.io import fits
from time import strftime
from tabulate import tabulate

""" make composite spectral type template for each subtype using
    spectra (R~10,000) from the ELODIE 3.1 stellar library.
    For each subtype, all existing spectra will be splined and coadded
    to create a composite template. When more than one spectrum exists
    for a star, it will be coadded and the resultant spectrum will be
    used to create the composite.
    
    Optionally, spectra of +/-1 subtype could be used to create the
    composite template.
    
    Written: October 2013 (SD)
"""


""" FIX: when combining, use median over mean
    : when N< 5 for a spT, using +-1 is fine
"""
    
def main():
    ## read in the elodie metadata
    elodie_path = '/Users/saurav/tara/local/datalib/spectral_templates/elodie/'
    elodie = np.genfromtxt(elodie_path+'table_meas.dat',names=True,usecols=range(10),
                           dtype='S8,S8,S8,f8,i2,f4,f8,i2,f4,i2,f4,i2')
    elodie = np.sort(elodie,order=['spT','name','num','Vr'])
    
    verbose = 2 # level of comments printed out: 0 --> none, 1 --> basics, 2 --> all, 3 --> all + plot
    
    threshdRV = 5. # threshold for dRV when a star has >1 spectra in   library
    spT = ['A0V',
           'F0V', 'F1V', 'F2V', 'F3V', 'F4V', 'F5V', 'F6V', 'F7V', 'F8V', 'F9V',
           'G0V', 'G1V', 'G2V', 'G3V', 'G4V', 'G5V', 'G6V', 'G7V', 'G8V', 'G9V',
           'K0V', 'K1V', 'K2V', 'K3V', 'K4V', 'K5V', 'K7V',
           'G0IV','G5IV','G8IV','G9IV','K1IV','K3IV',
           'G8III','K0III','K1III','K2III','K4III','K5III',
           'G0II']
#     spT = ['K3IV']
    nSpT = np.size(spT) 

    temSpec = None
    for i in range(nSpT):                
        # find the elodie stars of that spT
        subElodie = elodie[np.where((elodie['spT'] == spT[i]) &
                                    (elodie['name'] != 'SUN'))]
        nSpec = np.size(subElodie)
        
        if nSpec == 0:
            if verbose >= 1:
                print 'Warning: No match found in Elodie for %s. No template written.'%(spT[i])
        else:
            uniqNames, indx = np.unique(subElodie['name'],return_index=True)
            nStar = np.size(indx)
            
            if verbose >= 1:
                print '%s: No. of matching stars, spectra found: %i, %i' %(spT[i],nStar,nSpec)
            if verbose >= 2: # pretty print the elodie matches
                print tabulate(subElodie,headers=subElodie.dtype.names)
                print ''
            
        # find the unique stars of that spT
#         if nSpec < 10:
#             moreSpec = findMoreSpec(spT[i])
            objSpec = None
            # now for each unique star, loop
            for j in range(nStar):
                ind = np.where(subElodie['name'] == uniqNames[j])
                nMult = np.size(ind)
                      
                dRV = np.ptp(subElodie['Vr'][ind])
                if dRV > threshdRV:
                    if verbose >= 1:
                        print '\t%s: %s ignored. delta_RV > %.2f km/s.' %(subElodie['spT'][ind][0],
                                                                    subElodie['name'][ind][0],
                                                                    threshdRV)
                else:
                    # read and coadd (if nMult > 1) spectra 
                    for k in range(nMult):
                        spec = openSpec(elodie_path+'H/'+subElodie['num'][ind][k]+'.fits')
                      
                        # add the multiple spectra for that object -- equal weighting
                        objSpec = spec if objSpec == None else (objSpec + spec)
                
                    objSpec = normalize(objSpec)
                    
                    if verbose >= 3:
                        plt.xlabel('Wavelength')
                        plt.ylabel('Flux')
                        plt.axis([3800,7000,0,nStar+2])
                        plt.title(spT[i])
                    
                        plt.plot(objSpec.wave,(j+1)+objSpec.flux,color='blue')
                            
            # now, add the spectra for each star of that spT -- equal weighting
            temSpec = objSpec if temSpec == None else (temSpec + objSpec)
    
            # plot, and save spectrum
            if verbose >= 3:
                temSpec = normalize(temSpec) # normalize to 1
                plt.plot(temSpec.wave,temSpec.flux,color='red')
                fout = '%s.pdf' %(spT[i])
                plt.savefig(fout)
            
            # create a hdu to save temSpec
            hdu = fits.PrimaryHDU([spec.wave,spec.flux])
            # update headers
            hdu.header['CRVAL1'] = spec.wave[0]
            hdu.header['CD1_1'] = spec.wave[1] - spec.wave[0]
            hdu.header['spT'] = (spT[i],'Spectra Type of template')
            hdu.header['nSpectra'] = ('nSpec','No. of templates used to create composite')
            hdu.header['nStars'] = ('nStar','No. of unique stars used to create composite')
            hdu.header['Comment'] = ('Composite template created by coadding spectra Elodie3.1 Stellar Library')
            hdu.header['History'] = ('Created on: %s') %strftime("%Y-%m-%dT%H:%M:%S")
            # write the hdu to file
            hdulist = fits.HDUList([hdu])
            hdulist.writeto(spT[i]+'.fits', clobber=True, checksum=True)
            hdulist.close()
            
# read in spectrum; returns spec
def openSpec(filename):
    hdu = fits.open(filename)
    flux = hdu[0].data
    hdr = hdu[0].header
    wave = hdr['CRVAL1'] + hdr['CD1_1'] * np.arange(hdr['NAXIS1'])
    hdu.close()
    
    spec = psyn.ArraySpectrum(wave=wave,flux=flux)
    return spec

def normalize(spec):
    f = spec.flux
    w = spec.wave
    
    indx = np.where((w > 4500) & (w < 4550))
    f /= np.mean(f[indx])
    return psyn.ArraySpectrum(wave=w,flux=f)

def findMoreSpec(sp):
    
    return fname

def resample(spec,newWave):
    return 'resample me if you want to but elodie is on the same wavelength scale'
    
main()    
