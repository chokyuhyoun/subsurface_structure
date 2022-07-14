# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:34:21 2020

@author: chokh
"""

import os
import numpy as np
import copy
from glob import glob
from scipy import fftpack
import matplotlib.pyplot as plt
from mylib import misc
from astropy.io import fits
from astropy.time import Time
import mylib.image_process as ip


def filtering(data, dt, i_min=1, f_min=8):
    shape = data.shape
    freq_high = 1/60/i_min
    freq_low = 1/60/f_min
    fft_res = fftpack.fftn(data)
    freq_z = fftpack.fftfreq(shape[0], dt)[:, None, None]
    fil_band = (abs(freq_z) > freq_low) & (abs(freq_z) < freq_high)
    ndata = (fftpack.ifftn(fft_res*fil_band)).real
    return ndata

# y_init, x_init
init = np.array([[100, 100], \
                 [100, 170], \
                 [120, 120], \
                 [200, 340], \
                 [110,  60], \
                 [110,  80], \
                 [200, 345], \
                 [150, 350], \
                 [  0,   0], \
                 [125, -25]])
mar = np.array([0, 25, 0, 50, 0, 0, 0, 0])
mar_t = np.array([[ 0,  -1], \
                  [60,  -1], \
                  [60,  -1], \
                  [60,  -1], \
                  [ 0,  -1], \
                  [ 0,  -1], \
                  [50,  -1], \
                  [60,  -1], \
                  [ 0, 100], \
                  [ 0, 160]])
wid = 200                     
path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)

files_2796 = glob(path+'\iris*SJI_2796*.fits')
files_2832 = glob(path+'\iris*SJI_2832*.fits')

ref_images = np.zeros((len(files_2796), wid, wid))
for num in range(len(files_2796)):
    print(num)
    data_2796 = fits.getdata(files_2796[num])[mar_t[num, 0]:mar_t[num, 1]]
    hdr_2796 = fits.getheader(files_2796[num])
    data_2832 = fits.getdata(files_2832[num])[mar_t[num, 0]:mar_t[num, 1]]
    hdr_2832 = fits.getheader(files_2832[num])
    shape = np.array(data_2832.shape)
    cdelt_2796 = hdr_2796['CDELT1']
    dt_2796 = hdr_2796['CDELT3']
    nf = data_2796.shape[0]
    pdata_2796 = np.zeros([nf, wid, wid])
    pdata_2832 = np.zeros([nf, wid, wid])
    delyx = init[num]
    if num == 0:
        pre_2796 = data_2796[:, delyx[0]:delyx[0]+wid, 
                                delyx[1]:delyx[1]+wid] 
        for i in range(nf):
            pdata_2796[i] = ip.cr_mem(pre_2796[i])
        pdata_2832 = data_2832[:, delyx[0]:delyx[0]+wid, 
                                  delyx[1]:delyx[1]+wid] 
    else:        
        obj_image = data_2832[0, delyx[0]:delyx[0]+wid, 
                                 delyx[1]:delyx[1]+wid]
        delyx1, cor = ip.alignoffset(obj_image, ref_image, cor=True)
        delyx = delyx + delyx1.flatten().astype(int)
        
        for i in range(nf):
            pre_2796 = data_2796[i, delyx[0]:delyx[0]+wid, 
                                    delyx[1]:delyx[1]+wid]
            pdata_2796[i] = ip.cr_mem(pre_2796)
            pdata_2832[i] = data_2832[i, delyx[0]:delyx[0]+wid, 
                                         delyx[1]:delyx[1]+wid]
    ndata_2796 = filtering(pdata_2796, dt_2796)
    ref_image = pdata_2832[-1]

    filename_2796 = 'fp_2796_'+os.path.basename(files_2796[num])[0:23]+'.fits'
    filename_2832 = 'fp_2832_'+os.path.basename(files_2832[num])[0:23]+'.fits'
    nhdr_2796 = copy.deepcopy(hdr_2796)
    nhdr_2796['CRPIX1'] = 0.5*(wid-1)
    nhdr_2796['CRPIX2'] = 0.5*(wid-1)
    nhdr_2796['CRVAL1'] = hdr_2796['CRVAL1']- \
                            (hdr_2796['CRPIX1']-(delyx[1]+0.5*(wid-1)))*hdr_2796['CDELT1']
    nhdr_2796['CRVAL2'] = hdr_2796['CRVAL2']- \
                            (hdr_2796['CRPIX2']-(delyx[0]+0.5*(wid-1)))*hdr_2796['CDELT2']
    t0 = Time(hdr_2796['DATE_OBS']).jd
    nhdr_2796['DATE_OBS'] = Time(t0+hdr_2796['CDELT3']*mar_t[num, 0]/86400., format='jd').isot
    nhdr_2832 = copy.deepcopy(hdr_2832)
    nhdr_2832['CRPIX1'] = 0.5*(wid-1)
    nhdr_2832['CRPIX2'] = 0.5*(wid-1)
    nhdr_2832['CRVAL1'] = hdr_2832['CRVAL1']- \
                            (hdr_2832['CRPIX1']-(delyx[1]+0.5*(wid-1)))*hdr_2832['CDELT1']
    nhdr_2832['CRVAL2'] = hdr_2832['CRVAL2']- \
                            (hdr_2832['CRPIX2']-(delyx[0]+0.5*(wid-1)))*hdr_2832['CDELT2']
    t0 = Time(hdr_2832['DATE_OBS']).jd
    nhdr_2832['DATE_OBS'] = Time(t0+hdr_2832['CDELT3']*mar_t[num, 0]/86400., format='jd').isot

    fits.writeto(filename_2796, ndata_2796, nhdr_2796, overwrite=True)
    fits.writeto(filename_2832, pdata_2832, nhdr_2832, overwrite=True)
 
    ref_images[num] = copy.deepcopy(ref_image)
