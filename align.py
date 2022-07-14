# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:34:21 2020

@author: chokh
"""

import os
import numpy as np
import copy
# import astropy.units as u
import glob
# import fisspy.image.base as fp
from scipy import fftpack
from scipy import interpolate
import matplotlib.pyplot as plt
from mylib import misc
from astropy.io import fits
from astropy.time import Time, TimeDelta
import mylib.image_process as ip
import copy


# y_init, x_init
init = np.array([[100, 100], \
                 [100, 170], \
                 [120, 120], \
                 [200, 340], \
                 [110, 60], \
                 [110, 80], \
                 [200, 345], \
                 [150, 350]])
mar = np.array([0, 25, 0, 50, 0, 0, 0, 0])
mar_t = np.array([0, 60, 0, 60, 0, 0, 0, 0])
wid = 200                     
path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)

num = 0
files_2796 = glob.glob(path+'\iris*SJI_2796*.fits')
files_2832 = glob.glob(path+'\iris*SJI_2832*.fits')
data_2796 = fits.getdata(files_2796[num])[mar_t[num]:]
hdr_2796 = fits.getheader(files_2796[num])
data_2832 = fits.getdata(files_2832[num])[mar_t[num]:]
hdr_2832 = fits.getheader(files_2832[num])
shape = np.array(data_2796.shape)
cdelt_2796 = hdr_2796['CDELT1']
dt_2796 = TimeDelta(hdr_2796['cdelt3'], format='sec')
t_2796 = Time(hdr_2796['date_obs'])+dt_2796*np.arange(shape[0])
cdelt_2832 = hdr_2832['CDELT1']
dt_2832 = TimeDelta(hdr_2832['cdelt3'], format='sec')
t_2832 = Time(hdr_2832['date_obs'])+dt_2832*np.arange(shape[0])
t_2832_sec = (t_2832.jd - t_2832[0].jd)*86400
t_2796_sec = (t_2796.jd - t_2832[0].jd)*86400
nf = data_2796.shape[0]

pdata_2796 = np.zeros([nf, 200, 200])
pdata_2832 = np.zeros([nf, 200, 200])
del_2796 = np.zeros([nf, 2])
del_2832 = np.zeros([nf, 2])

yyp0, xxp0 = np.mgrid[0:shape[1], 0:shape[2]]
yyp, xxp = np.mgrid[0:wid, 0:wid]    

for i in range(nf):
    if (i == 0) & (num == 0):
        del_2796[i] = init[num]
        del_2832[i] = init[num]
        ref_2796 = data_2796[i, init[num, 0]:init[num, 0]+wid, 
                                init[num, 1]:init[num, 1]+wid]
        pdata_2796[i, :, :] = copy.deepcopy(ref_2796)
        ref_2832 = data_2832[i, init[num, 0]:init[num, 0]+wid, 
                                init[num, 1]:init[num, 1]+wid]
        pdata_2832[i, :, :] = copy.deepcopy(ref_2832)
        continue
    if (i == 0) & (num != 0):
        del_2796[i-1] = init[num]
        del_2832[i-1] = init[num]
        dum1 = glob.glob(path+r'\fp_2832*.fits')
        dum2 = fits.getdata(dum1[num-1])
        ref_2832 = copy.deepcopy(dum2[-1, :, :])

    pdata_2832t = ip.interpol2d(data_2832[i], xxp+del_2832[i-1, 1], yyp+del_2832[i-1, 0])

    # ref_2832[ref_2832 > 700] = 700
    # ref_2832[ref_2832 < 0] = 0
    # pdata_2832t[pdata_2832t > 700] = 700
    # pdata_2832t[pdata_2832t < 0] = 0
    
    delyx, cor = ip.alignoffset(pdata_2832t[:, mar[num]:-mar[num]-1], 
                                ref_2832[:, mar[num]:-mar[num]-1], cor=True)
    if (cor >= 0.5) | (i == 0):
        delyx = delyx.flatten()
        del_2832[i] = del_2832[i-1] + delyx
        yyp1 = yyp + del_2832[i, 0]
        xxp1 = xxp + del_2832[i, 1]
        ref_2832 = ip.interpol2d(data_2832[i], xxp0, yyp0, xxp1, yyp1)
    else:
        print(i)
        ref_2832 = copy.deepcopy(pdata_2832[i-1])
        del_2832[i] = copy.deepcopy(del_2832[i-1])
    pdata_2832[i, :, :] = copy.deepcopy(ref_2832)

    if (i == 0) & (num != 0):
        ref_2796 = ip.interpol2d(data_2796[i], xxp0, yyp0, xxp1, yyp1)
        pdata_2796[i] = copy.deepcopy(ref_2796)
        del_2796[i] = copy.deepcopy(del_2832[i])
        continue
    
    pdata_2796t = ip.interpol2d(data_2796[i], xxp+del_2796[i-1, 1], 
                                              yyp+del_2796[i-1, 0])
    ref_2796[ref_2796 <= 0] = 0
    ref_2796[ref_2796 > 100] = 100
    pdata_2796t[pdata_2796t <= 0] = 0
    pdata_2796t[pdata_2796t > 100] = 100
    delyx, cor = ip.alignoffset(pdata_2796t[:, mar[num]:-mar[num]-1], 
                                ref_2796[:, mar[num]:-mar[num]-1], cor=True)
    del_2796[i] = del_2796[i-1] + delyx.flatten()
    
    yyp1 = yyp + del_2796[i, 0]
    xxp1 = xxp + del_2796[i, 1]
    pdata_2796[i] = ip.interpol2d(data_2796[i], xxp0, yyp0, xxp1, yyp1)
    ref_2796 = copy.deepcopy(ip.cr_mem(pdata_2796[i]))
    # pdata_2796[i, :, :] = ip.cr_mem(ref_2796, size=1)
    # print(i, del_2796[i], del_2832[i])
dt = hdr_2796.get('cdelt3')  # in second
freq_high = 1/60/1
freq_low = 1/60/8
t = np.arange(0, nf*dt, dt)
fft_res = fftpack.fftn(pdata_2796)

freq_z = fftpack.fftfreq(nf, dt)[:, None, None]
fil_band = (abs(freq_z) > freq_low) & (abs(freq_z) < freq_high)
ndata_2796 = (fftpack.ifftn(fft_res*fil_band)).real

filename_2796 = 'fp_2796_'+os.path.basename(files_2796[num])[0:23]+'.fits'
filename_2832 = 'fp_2832_'+os.path.basename(files_2832[num])[0:23]+'.fits'
nhdr_2796 = copy.deepcopy(hdr_2796)
nhdr_2796['CRPIX1'] = hdr_2796['CRPIX1']-init[num, 1]
nhdr_2796['CRPIX2'] = hdr_2796['CRPIX2']-init[num, 0]
fits.writeto(filename_2796, ndata_2796, nhdr_2796, overwrite=True)
fits.writeto(filename_2832, pdata_2832, hdr_2832, overwrite=True)


