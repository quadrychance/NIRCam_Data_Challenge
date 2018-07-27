from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sep
from astropy.io import ascii
from astropy import wcs
from astropy.stats import sigma_clipped_stats 
from matplotlib.patches import Ellipse



mask = np.genfromtxt('positions.txt')



data = fits.getdata('test_mosaic_F200W_2018_03_19.fits', ext=0)
data = data.byteswap().newbyteorder()
bkg = sep.Background(data, bw = 32)
back = bkg.back() 
rms = bkg.rms()
data_sub = data - bkg


thresh = 20 * bkg.globalrms
objects = sep.extract(data_sub, thresh)



bkg_image = bkg.back()

plt.imshow(bkg_image, interpolation='nearest', origin='lower')
#plt.colorbar();
plt.show()
bkg_rms = bkg.rms()
plt.imshow(bkg_rms, interpolation='nearest',  origin='lower')
#plt.colorbar();
plt.show()
n = np.asarray(bkg_rms)
hdu = fits.PrimaryHDU(n)
hdu.writeto('rms.fits')

coords = (objects['x'], objects['y'])
coords2 = np.array( (objects['x'], objects['y']))
coordsT = coords2.T
print(np.shape(coords))
np.savetxt('objects.txt', coords, fmt='%.3f', delimiter = '	')


fig, ax = plt.subplots()
m, s = np.mean(data_sub), np.std(data_sub)
im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')


for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=6*objects['a'][i],
                height=6*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
theta=0
b =1
apertures = [2.5,3.4,5,7.5,10,15,20,30,50,70]
catalog = []
f = []
a =1
flux, fluxerr, flag = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[0],err=bkg.globalrms, gain=1.4)
flux1, fluxerr1, flag1 = sep.sum_ellipse(data_sub, objects['x'], objects['y'],  a, b, theta,apertures[1],err=bkg.globalrms, gain=1.4)
flux2, fluxerr2, flag2 = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[2],err=bkg.globalrms, gain=1.4)
flux3, fluxerr3, flag3 = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[3],err=bkg.globalrms, gain=1.4)
flux4, fluxerr4, flag4 = sep.sum_ellipse(data_sub, objects['x'], objects['y'],  a, b, theta,apertures[4],err=bkg.globalrms, gain=1.4)
flux5, fluxer5r, flag5 = sep.sum_ellipse(data_sub, objects['x'], objects['y'],  a, b, theta,apertures[5],err=bkg.globalrms, gain=1.4)
flux6, fluxerr6, flag6 = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[6],err=bkg.globalrms, gain=1.4)
flux7, fluxerr7, flag7 = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[7],err=bkg.globalrms, gain=1.4)
flux8, fluxerr8, flag8 = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[8],err=bkg.globalrms, gain=1.4)
flux9, fluxerr9, flag9 = sep.sum_ellipse(data_sub, objects['x'], objects['y'], a, b, theta,apertures[9],err=bkg.globalrms, gain=1.4)
r, flag = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r1, flag1 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r2, flag2 = sep.flux_radius(data_sub, objects['x'], objects['y'], 6.*objects['a'],0.5, normflux=flux, subpix=5)
r3, flag3 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r4, flag4 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r5, flag5 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r6, flag6 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r7, flag7 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r8, flag8 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
r9, flag9 = sep.flux_radius(data_sub, objects['x'], objects['y'],6.*objects['a'], 0.5, normflux=flux, subpix=5)
    
#print(flux)
    #return flux
          

catalog = np.array(flux)
catalog1 = np.array(flux1)
catalog2 = np.array(flux2)
catalog3 = np.array(flux3)
catalog4 = np.array(flux4)
catalog5 = np.array(flux5)
catalog6 = np.array(flux6)
catalog7 = np.array(flux7)
catalog8 = np.array(flux8)
catalog9 = np.array(flux9)



half_light = np.array(r)
half_light1 = np.array(r1)
half_light2 = np.array(r2)
half_light3 = np.array(r3)
half_light4 = np.array(r4)
half_light5 = np.array(r5)
half_light6 = np.array(r6)
half_light7 = np.array(r7)
half_light8 = np.array(r8)
half_light9 = np.array(r9)



hdulist = fits.open('test_mosaic_F200W_2018_03_19.fits')
w = wcs.WCS(hdulist[0].header)
world = w.wcs_pix2world(coordsT,0)
worldT = world.T



master_catalog = (np.vstack((coords,worldT,catalog, catalog1, catalog2, catalog3, catalog4,  catalog5 ,catalog6, catalog7, catalog8,catalog9, half_light8))).T
#ascii.write(master_catalog, 'master_catalog.txt', overwrite = 'True', format = 'tab') 
print(np.shape(master_catalog))


np.savetxt('master_catalog', master_catalog, fmt='%.2f', delimiter = '	', header = 'x	y	RA	DEC	2.5	3.4	5.0	7.5	10	15	20	30	50	70	half light r')
plt.show()
