""" Perform aperture photometry on the reference of images given x,y coords

    Make sure the input, telescope, and photomotery params are correct.
    They are all specified at the beginning of this file.

    This is a special case of aper.py and plots the given stars (x, y) for 
    verification of image registration. Use this for testing purposes.

    Written by: S Dhital (June 2015)
"""
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from matplotlib.colors import LogNorm
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.table import hstack
from tabulate import tabulate

baseFname = 'ref.' + 'PG1613+426293'
imFile =  baseFname + '.fits'  # image file name
coordFile = baseFname + '.coords' # list of coordinates (x, y)
photFile = baseFname + '.phot' # output filename

aper_size = 8.
annulus, dannulus = 15., 20.

# read the reference image
im = pyfits.getdata(imFile, ignore_missing_end=1)

# read in coordinates of the stars -- assumes (N,2)
coords = np.genfromtxt(coordFile, dtype='f8')

print coords.shape

# max, min pixel values to be plotted
vmin, vmax = 0.1, 70000
# normalization to be used -- add more options
norm = LogNorm(vmin, vmax)
# norm=ImageNormalize(stretch=SqrtStretch())

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.imshow(im, vmin=vmin,vmax=vmax, cmap='Greys', norm=norm)

for i in range(np.size(coords[:,0])):
    circ = plt.Circle((coords[i,1:2]), aper_size, color='r',fill=False)
    ax.add_patch(circ)

#plt.savefig(baseFname + '.png',format='png')
plt.show()

""" Aperture photometry happens here.
    follow the algorithm outlined in 
    http://photutils.readthedocs.org/en/latest/photutils/aperture.html#global-background-subtraction
    but smooth background error not included for now
    """
# create aperture and annulus objects
apertures = CircularAperture(coords, r=aper_size)
annulus_apertures = CircularAnnulus(coords, r_in=annulus, r_out=dannulus)
   
# calculate the object and annulus flux
data_error = np.sqrt(im)
rawflux_table = aperture_photometry(im, apertures, error=data_error)
bkgflux_table = aperture_photometry(im, annulus_apertures, error=data_error)
# stack the two tables into one
phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])

# calculate bkg mean -- normalize by area
bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
bkg_sum = bkg_mean * apertures.area()

# final photometric flux
final_sum = phot_table['aperture_sum_raw'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum

#print phot_table.columns
np.savetxt(photFile, phot_table, fmt='%10.2f' )
