import numpy as np
from astropy.io import ascii
from astropy.io import fits
import matplotlib.pyplot as plt
from pylab import rcParams

#script for obtaining spectra from the JADES Mock Catalog simulated images
#written in python 3.6
#to run enter python mock_catalog_spectra.py into the terminal



rcParams['figure.figsize'] = 15, 15

np.set_printoptions(threshold='nan')

#open fits files

hdul2 = fits.open('JADES_SF_mock_r1_v1.1_spec_5A_ID_200001_250000.fits')
hdul3 = fits.open('JADES_SF_mock_r1_v1.1_spec_5A_ID_250001_300000.fits')


#crack the fits files into the sets of fits metadata that contains the fluxes and the wavelength range
#extension 1 contains the fluxes for JADES Mock Catalog files
#extensions 2 and 3 contain the wavelength and ID/ redshifts data respectively but we only need the wavelength
spec = hdul[1].data
spec2 = hdul2[1].data
spec3 = hdul3[1].data

wavelength = hdul[2].data

#since the data is broken up over multiple fits files, the spectra of any particular galaxy needs to be matched with its number in the array containing them.
#For example, the spectra of object 188001 is 38001 entries deep into the numpy array representing the fits extension beginning with object 150001
spec_188001 = spec[38001,:]
spec_175101 = spec[25101,:]
spec_189701 = spec[39701,:]




spec_217901 = spec2[17901,:]
spec_211201 = spec2[11201,:]
spec_203901 = spec2[3901,:]
spec_206001 = spec2[6001,:]


spec_263701 = spec3[13701,:]
spec_284101 = spec3[34101,:]
spec_271801 = spec3[21801,:]
spec_288701 = spec3[38701,:]
spec_271401 = spec3[21401,:]






#the above process yields a bunch of 1-d arrays that need to be concatenated along the short axis
spec_com = np.vstack((spec_188001,spec_175101,spec_189701,spec_217901,spec_211201,spec_203901,spec_206001,spec_263701,spec_284101,spec_271801,spec_288701,spec_271401))

#for whatever reason this array needs to be transposed
spec_comT = spec_com.T

#making sure the wavelength array is of the proper dimensions
wav = np.reshape(wavelength,(3075,1))

#combine wavelentgh and fluxes into 2d array that is easily plottable
spec_data = np.hstack((wav,spec_comT))#spec_188001,spec_175101,spec_189701,spec_217901,spec_211201,spec_203901,spec_206001,spec_263701,spec_284101,spec_271801,spec_288701,spec_271401) )

#the rest is just plotting.
#plt makes the subplots interact weirdly and isn't as amenable to being put in a loop so ax objects are used here
#
fig, ax = plt.subplots(12, 1, sharex='col')
fig.subplots_adjust(hspace=1)

ax[6].set_ylabel('$erg cm^{-2} s^{-1}}$', size = 'x-large')
ax[11].set_xlabel('$\AA$', size = 'x-large' )

names = ['1','spec_188001','spec_175101','spec_189701','spec_217901','spec_211201','spec_203901','spec_206001','spec_263701','spec_284101','spec_271801','spec_288701','spec_271401']

for i in range(1,12):
    ax[i].plot(spec_data[:,0],spec_data[:,i])
    #ax[i].text(0.5, 0.5, names[i], fontsize=18, ha='center')
    ax[i].set_title(names[i])
    ax[i].set_yscale('log')
    ax[i].set_xscale('log')
   


plt.savefig('spectra', dpi = 1000)
#!!!!!LOOK HERE!!!!!!
# for whatever reason, the formatting of the spec_comT array is NOT preserved when saving using np.savetxt or getting it directly from stdout
ascii.write(spec_data,'specdata.txt', overwrite = 'True',delimiter = ' ', names=('wavelength', 'spec_188001','spec_175101','spec_189701','spec_217901','spec_211201','spec_203901','spec_206001','spec_263701','spec_284101','spec_271801','spec_288701','spec_271401'),formats ={'wavelength':'%.4g', 'spec_188001':'%.4g','spec_175101':'%.4g','spec_189701':'%.4g','spec_217901':'%.4g','spec_211201':'%.4g','spec_203901':'%.4g','spec_206001':'%.4g','spec_263701':'%.4g','spec_284101':'%.4g','spec_271801':'%.4g','spec_288701':'%.4g','spec_271401':'%.4g'})















