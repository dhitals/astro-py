""" Calculates the Scargle periodogram as per the algorithm in Horne & Baliunas (1986)

    INPUT:  dstar###.f.dat
    OUTPUT: pstar###.f.day

    HISTORY:- Written by S. Dhital -- July 2014

    FIX: j does not reset when moving to the next object
"""

import numpy as np
import os

nIter = 300000 # no. of steps in the Scargle periodogram
omega = 2 * np.pi * 0.005 # initial period = 0.005 days
    
starList = [file for file in os.listdir('.') 
    if file.endswith ('dstar',0,5) if file.endswith('.dat')]
print starList, np.size(starList)

## do for each star
for i in range(np.size(starList)):
    #read in the dstar###.dat file: expects: ID, jdate, dmag, magerr
    star = np.genfromtxt(starList[i], names=True, dtype='S3, f8, f8, f8')

    nObs = np.size(star)
    # rangeJD = max(star['jdate']) - min(star['jdate'])
    avgMag = sum(star['dmag']) / nObs
    var = sum( (star['dmag'] - avgMag)**2 / (nObs-1) )  # what is the difference between stddev and this?
    
    # make an array to store all the output info
    dtype = [('period','f8'),('power','f8'),('fap','f8')]
    outArray = np.zeros(nIter,dtype=dtype)
        
    for j in range(nIter):
        s5 = sum( np.sin(2*omega*star['jdate']) )
        s6 = sum( np.cos(2*omega*star['jdate']) )
        tau = (np.arctan2(s5,s6) / (2*omega)) # tau from Eq (2) in H&B86
            
        s1 = sum( (star['dmag'] - avgMag) * np.cos(omega * (star['jdate']-tau) ) )
        s2 = sum( np.cos( omega*(star['jdate'] - tau) )**2 )
        s3 = sum( (star['dmag'] - avgMag) * np.sin(omega * (star['jdate']-tau) ) )
        s4 = sum( np.sin( omega * (star['jdate']-tau) )**2 )
            
        Px = 0.5 * ( (s1**2/s2) + (s3**2/s4) )  # Eq (1) in H&B86
        power = Px/var
            
        Ni = -6.362 + 1.193*nObs + 0.00098*nObs**2 # H&B fap -- likely not useful
        fap = 1 - (1 - np.exp(-1*power))**Ni
            
        period = omega / (2 * np.pi)
            
        outArray['period'][j] = period
        outArray['power'][j] = power
        outArray['fap'][j] = fap

        ## imcrease omega based on current value
        omega += (5e-4 if period > 1.5 else 5e-3)
        
        fout = 'p'+starList[i][1:] # pstar###.f.dat
        hdr = 'period    power    fap'
        fmt = '%1.5f'
        np.savetxt(fout, outArray, fmt=fmt, header=hdr)
