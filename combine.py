# -*- coding: utf-8 -*-
"""
Quadry Chance
Python 3.6

Combining JADES and GALAPAGOS catalogs
"""

import numpy as np
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.units import cds
from astropy.table import Table, join, hstack
from astropy.coordinates import match_coordinates_sky, concatenate, search_around_sky
import matplotlib.pyplot as plt
from astropy.io import fits

#import results from catstrip.py
hdul = fits.open('JADES_SF_mock_r1_v1.1.fits')
jadescat = hdul[1].data
hdul2 = fits.open('combcattest.fits') #fits file output from galapagos
galcat = hdul2[1].data
print("reading catalog files...")
ascii.write(galcat,'galcat.txt', overwrite = 'True',delimiter = '	')
ascii.write(jadescat,'test.txt', overwrite = 'True',delimiter = ' ')
galcat = np.genfromtxt('galcat.txt', skip_header = 1)
jadescat = np.genfromtxt('jadescat.txt', skip_header=1)


#taking only the needed columns from each
gal_short = np.vstack((galcat[:,0], galcat[:,12], galcat[:,13], galcat[:,42],galcat[:,43], galcat[:,40]*0.0317, galcat[:,41]*0.0317, galcat[:,38], galcat[:,39]))
jades_short = np.vstack((jadescat[:,0],jadescat[:,1], jadescat[:,2],jadescat[:,78],jadescat[:,75], (-2.5)*(np.log10(jadescat[:,24])+8.9)))
#truncate the floats within
gal2 = np.round(gal_short,5)
jades2 =np.round(jades_short,5)

#convert numpy arrays to astropt.Table objects
print('converting catalog to table...')
g = Table(gal2.T, names=('CATALOG_#','RA','DEC', 'GALFIT_SERSIC_N','GALFIT_SERSIC_NERR', 'GALFIT_R_E','GALFIT_R_EERR', 'GALFIT_MAG', 'GALFIT_MAGERR'))
j = Table(jades2.T, names=('CATALOG_#','RA','DEC', 'SERSIC_N', 'R_E', 'MAG'))



#convert RA and DEC in both catalogs to SkyCoord objects in order to use astropy.coordinates.matching.search_around_sky  
print('converting to SkyCoord...')
gal_coord = SkyCoord(ra=g['RA']*u.degree, dec=g['DEC']*u.degree)  
jades_coord = SkyCoord(ra=j['RA']*u.degree, dec=j['DEC']*u.degree) 


max_sep = .1*cds.arcs #maxiumum separation of matched galaxies

#match galaxies in galcat to galaxies within jades cat. 
#only match the smaller catalog into the larger to avoid multiple matches per object
#idx1 and idx2 are indexes within jades and gal, respectively
print('matching objects...')
idx1, idx2, d2d, d3d = gal_coord.search_around_sky(jades_coord, seplimit = max_sep) 
separation = np.asarray(d2d) #separation from nearesr source in the other catalog

#index catalogs by matches and concatenate separations onto the array
jades_matches = np.column_stack((jades_short.T[idx1], separation.T))
gal_matches = np.column_stack((gal_short.T[idx2], separation.T))


#make arrays into tables and combine them
print('putting matches into table')
jades_table  = Table(jades_matches, names=('CATALOG_#','RA','DEC', 'SERSIC_N', 'R_E', 'MAG', 'SEPARATION'))
gal_table  = Table(gal_matches, names=('CATALOG_#2','RA','DEC', 'GALFIT_SERSIC_N','GALFIT_SERSIC_NERR', 'GALFIT_R_E','GALFIT_R_EERR', 'GALFIT_MAG', 'GALFIT_MAGERR', 'SEPARATION'))                               
comb_table = hstack([jades_table,gal_table])


#write table to file
#for a program readable table, use format = 'tab'                            
ascii.write(comb_table,'comb_table.txt', overwrite = 'True', format = 'fixed_width')

'''plt.hist(jades_table['SEPARATION'], bins = 10)
plt.hist2d(jades_table['RA'],jades_table['SEPARATION'], bins = 100, cmap = 'inferno')
plt.hist2d(jades_table['DEC'],jades_table['SEPARATION'], bins = 100,cmap = 'inferno')

plt.hist(gal_table['SEPARATION'], bins = 100)
plt.hist2d(gal_table['RA'],gal_table['SEPARATION']*3600, bins = 100, cmap = 'inferno')
plt.hist2d(gal_table['DEC'],gal_table['SEPARATION'], bins = 100,cmap = 'inferno')

plt.ylim(0,0.0001)
plt.show()'''
print('done')
