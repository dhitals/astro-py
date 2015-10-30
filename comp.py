import numpy as np
from tabulate import tabulate
import os

"""
 INPUT:  phot.dat, dstar.lis, cstar.lis
 OUTPUT: dstar###.f.dat, sigma.f.out, cstarsigma.f.out

 HISTORY:- Written in FORTRAN by William Herbst et al.
           Written in IDL by Saurav Dhital  -- May 2006
           Written in python by Saurav Dhital -- July 2014

 Fix: when only one dstar, it crashes
    : probably better to not include own mag when comparing sigInt, sigExt
    : possibly print out avgmag when printing out resulting paramters
    : Modularize!
"""

### Check these parameters before running code

path = 'PG1613+426/reducedimages/' #where your images are at 
photFile = 'phot.dat'   # default photometry input file, where the phot file is at

filt = 'None' # Which filter are you working on?

magThresh = 25.0 # magnitude threshold: fainter mags will be ignored

# 0.050 is a good starting value, eventually you want to end up with < 0.015
magErrThresh = 0.10 # magErr threshold: errors larger than this value will be ignored

# print dstar files? Do this only after you are satisfied with your comp stars
print_dstar_file = 'yes'

### now read input files
# 'star' is a 2D array that stores all the input photometry info 
# expected cols: imageName, jdate, airmass, fwhm, starNum, x, y, mag, magerr
star = np.genfromtxt(photFile,names=True,dtype='S22,f8,f8,f8,S3,f8,f8,f8,f8')

# read in list of stars for which dmag is to be calculated -- can be the list of all stars
dstarList = np.genfromtxt('dstar.lis',dtype='S3',comments='#')
dnum = np.size(dstarList)

# read in the file with the list of comp stars -- initially the same as dstar
cstarList = np.genfromtxt('cstar.lis',dtype='S3',comments='#')
cnum = np.size(cstarList)    

# if only 1 comp star then use high magErrThreshold
if cnum == 1:
    magErrThresh = 9.0

# determine the no. of observations for the field -- based on unique JD
obsJD = np.unique(star['jdate'])
Nobs_tot = np.size(obsJD)

print ' '
print 'Filter selected: ',filt
print 'No. of unique epochs for this field: ',Nobs_tot
print 'No. of stars to be processed: ',dnum
print 'No. of comparison stars: ',cnum
print 'Magnitude threshold: ',magThresh
print 'Magnitude Error threshold: ',magErrThresh
print ' '    


""" calculations begin here  -- 
    this section calculates the mean mag for the comp stars, 
    so only works with a subset of the stars read in from photFile
    calculate <mag> only for (mag < magThresh && magerr < magErrThresh)
"""
# make a subset of cStars
cstar = star[np.in1d(star['starNum'],cstarList)] # select the subset of cstars
obsJD = np.unique(cstar['jdate'])
Nobs_cstar = np.size(obsJD)
print 'No. of unique epochs in which comp stars were observed: ', Nobs_cstar
if Nobs_cstar < Nobs_tot:
    print('\nWarning: Your comp stars were only observed in %i of %i epochs.\n' %(Nobs_cstar, Nobs_tot))

cMag = np.zeros(Nobs_cstar, dtype=[('jdate','f8'),('meanMag','f8'),('nStar','i8')])

# for each unique JD calculate the mean mags of comp stars
for i in range(Nobs_cstar):
    # select the ones for that JD & good mag, magerr
    indx = np.where((cstar['jdate'] == obsJD[i]) &
                    (cstar['mag'] < magThresh) &
                    (cstar['magerr'] < min(magErrThresh, 9)))
    
    cMag['jdate'][i] = obsJD[i] # store obsJD for this cmag for book-keeping

    if np.size(indx) > 0:    
        cMag['meanMag'][i] = sum(cstar['mag'][indx]) / np.size(indx)
        cMag['nStar'][i] = np.size(indx)
    else:
        noCstar = np.where(star['jdate'] == obsJD[i])
        #print('No comp stars w/ good mag/magErr found for %s. Removed from phot DB.' %(star['imageName'][noCstar][0]))

        star = star[np.where(star['jdate'] != obsJD[i])]
        #os.system( "mv %s %sbadshift/" %(path+s1['imageName'][0], path) )

# print histogram of how many times cstars were observed
hist = np.histogram(cMag['nStar'],bins=[0,1,2,3,4,5,6])
hist[0][0] = Nobs_tot - Nobs_cstar # add in the Ncount for N_images = 0

print tabulate((hist[1],hist[0]))
# pretty print JD, cmag, and num (num of comp stars for that image)
#print 'Here are mean mags of your comp stars for each epoch:'
#print tabulate(cMag,headers=['jdate','meanMag','nStar'],floatfmt='1.3f')

"""
Here we calculate the differential magnitude (dmag) for 
all the stars of interest w.r.t. your comp stars

unlike for cstar, here we iterate over each star (not each epoch)
"""
# make an array 'sigArray' to store all the output info
dtype = [('starNum','S3'),('NobsStar','i1'),('avgDmag','f8'),('rangeDmag','f8'),
         ('sigInt','f8'),('sigExt','f8'),('sigRatio','f8')]
sigArray = np.zeros(dnum, dtype=dtype)

# loop for eacgh dstar
for j in range(dnum):
    indx = np.where((star['starNum'] == dstarList[j]) &
                    (star['mag'] < magThresh) &
                    (star['magerr'] < min(magErrThresh, 9)) )
                   # (star['filter'] == filt))
    subStar = star[indx]
    Nobs_dstar = np.size(subStar) # no. of good images for that star

    # calculate the total dmag for each star
    dmag = [] #np.ones(Nobs_dstar)*99 # 99 because 0,1 are valid values for dmag
    for i in range(Nobs_dstar):
        indx2 = np.where(subStar['jdate'][i] == cMag['jdate']) # find the right epoch
        dmag.append(subStar['mag'][i] - cMag['meanMag'][indx2])

    # calculate the sigma info for each dstar -- for reference later
    sigArray['starNum'][j] = dstarList[j]
    sigArray['NobsStar'][j] = Nobs_dstar
    sigArray['avgDmag'][j] = avgDmag = sum(dmag) / Nobs_dstar
    sigArray['rangeDmag'][j] = max(dmag) - min(dmag)

    sigArray['sigInt'][j] = max( sum(subStar['magerr'])/Nobs_dstar, 0.007 ) # set a floor of 0.007
    sigArray['sigExt'][j] = np.sqrt( sum((dmag - avgDmag)**2) / (Nobs_dstar-1) )
    sigArray['sigRatio'][j] = sigArray['sigExt'][j] / sigArray['sigInt'][j]

    if Nobs_dstar < (0.75*Nobs):
        print('Warning: Star %s found in < 75%% (%i/%i) of your images.' %(dstarList[j], Nobs_dstar, Nobs))
    
    # only to be printed at the end -- when suitable comp stars have been found
    if print_dstar_file == 'yes':
        
        dtype = [('starNum', 'S3'), ('jdate','f8'), ('dmag','f8'), ('magerr','f8')]
        dmagFile = np.zeros(Nobs_dstar, dtype=dtype)
        
        print subStar.shape, dmagFile.shape, np.shape(dmag)

        dmagFile['starNum'] = subStar['starNum']
        dmagFile['jdate'] = subStar['jdate']
        dmagFile['dmag'] = dmag
        dmagFile['magerr'] = subStar['magerr']
        
        fout = 'dstar%3s.dat' %(dstarList[j])

        hdr = 'ID jdate          dmag  magerr'
        fmt = '%s %0.6f %8.3f %.3f'
        np.savetxt(fout, dmagFile, fmt=fmt, header=hdr)
            
        if dmagFile['starNum'][0] == '010': 
            print tabulate(dmagFile,headers=dmagFile.dtype.names)

## print out the results

# sigma_dstar.out  -- for all stars
fout = 'sigma_dstar.out'
cmt = '\n\t Sigma information for all stars (N = %s):\n' %(dnum)
hdr_file = 'ID  Nobs   dmag   rangeDmag  sigInt  sigExt  sigRatio'
fmt = '%3s %3d %8.3f %8.3f %8.3f %8.3f %8.3f'

# save to file
np.savetxt(fout,sigArray,fmt=fmt,header=hdr_file,comments=cmt,delimiter='\t')

# pretty-print to screen
hdr_scr = ['ID','Nobs','dmag','rangeDmag','sigInt','sigExt','sigRatio']
print cmt
print tabulate(sigArray, floatfmt='1.3f', headers=hdr_scr)

# sigma_cstar.out --- iff N_cstars > 1
if cnum > 1:
    indx3 = np.in1d(sigArray['starNum'][:],cstarList) # select the subset of cstars
    cstarSigArray = sigArray[indx3]
    
    fout = 'sigma_cstar.out'
    cmt = '\t Sigma information for comparison stars (N = %s):\n' %(cnum)
    ftr = 'Average sigExt for comparison stars is', sum(cstarSigArray['sigExt'][indx3])/cnum

    # sort the array for visual ease
    cstarSigArray = np.sort(cstarSigArray, order='sigInt')
    
    # save to file
    np.savetxt(fout,cstarSigArray,fmt=fmt,header=hdr_file,comments=cmt,delimiter='\t')

    # pretty print comp star info to screen
    print '\n'+cmt
    hdr = ['ID','NobsStar','dmag','rangeDmag','sigInt','sigExt','sigRatio']
    print tabulate(np.transpose(cstarSigArray),headers=hdr_scr,floatfmt='1.3f')
    print '\n \t %s %.3f'% (ftr)                            
