#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 14:06:58 2018

@author: quadrychance
creates rms image
"""

from astropy.io import fits
import numpy as np


data = fits.getdata('/home/quadrychance/galfit-example/EXAMPLE/test.fits', ext=0)
data = data.byteswap().newbyteorder()
rms = np.sqrt(abs(data))
n = np.array(rms)
hdu = fits.PrimaryHDU(n)
hdu.writeto('/home/quadrychance/galfit-example/EXAMPLE/rms.fits')
