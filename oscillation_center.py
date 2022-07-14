# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:36:27 2020

@author: chokh
"""

import os
import numpy as np
# import astropy.units as u
from glob import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.time import Time, TimeDelta
from mylib import nave, image_process, clump_find
from mylib import misc
import scipy.ndimage as spnd
import pickle


path = r'C:/Users/chokh/.spyder-py3/20151216'
os.chdir(path)
file_corr = glob(path+r'/clump_corr*.fits')
file_2796 = glob(path+r'/fp_2796*.fits')
file_2832 = glob(path+r'/fp_2832*.fits')
file_hmi = glob(path+r'/hmi*.fits')
hmi_data = fits.getdata(file_hmi[0])[:, 50:-50, 50:-50]
hmi_shape = np.array(hmi_data.shape)
with open(path+r'/hmi_time_info.p', 'rb') as file:
    hmi_time_str = pickle.load(file)
hmi_time = Time(hmi_time_str)

fig01, ax0n = plt.subplots(3, 1, figsize=[7, 9])
ax00 = fig01.add_axes([0, 0, 1, 1])
ax00.set_xlim(0, 1)
ax00.set_ylim(0, 1)
ax00.set_zorder(-1)
height = fig01.get_figheight()*fig01.get_dpi()
width = fig01.get_figwidth()*fig01.get_dpi()

ax01 = fig01.add_axes([0.12, 0.7, 0.3, 0.3])
ax01.set_xlabel('Solar X (arcsec)')
ax01.set_ylabel('Solar X (arcsec)')
ax01.set_xticks(np.arange(0, 150, 50))
ax01.set_yticks(np.arange(0, 150, 50))
hmi_extent = 0.5*np.array([0, hmi_shape[1], 0, hmi_shape[2]])
im01 = ax01.imshow(hmi_data[0], origin='lower', extent=hmi_extent, cmap='gray')
im01.set_clim(0, 6e4)
slit_ypos = 97
p01 = ax01.plot(hmi_extent[0:2], [slit_ypos*0.5]*2, 'r')
ax01.set_title('HMI Intensity')
t01 = ax01.text(10, 90, hmi_time[0].iso[0:-4]+' UT')
t01.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))

ax02 = fig01.add_axes([0.12, 0.1, 0.3, 0.55])
ax02.set_xlabel('Solar X (arcsec)')
ax02.set_ylabel('Time (UT)')
ax02.set_xticks(np.arange(0, 150, 50))
td_extent = np.array([hmi_extent[0:2], hmi_time.jd[[-1, 0]]]).flatten()
im02 = ax02.imshow(hmi_data[::-1, slit_ypos, :], origin='lower', extent=td_extent, 
                   cmap='gray')
im02.set_clim(im01.get_clim())
ax02.set_aspect('auto')
td_ytick_jd = np.fix(hmi_time[0].jd)+np.array([1.5, 2.5, 3.5])
td_yticklabel = ['2015 \n'+iso[5:10]+'\n00:00' \
                 for iso in Time(td_ytick_jd, format='jd').iso]
ax02.set_yticks(td_ytick_jd)
ax02.set_yticklabels(td_yticklabel)

for i, ax03 in enumerate(ax0n):
    num = i*2
    ax03.set_position([0.58, 0.675-i*0.32, 0.32, 0.32])
    ax03.set_aspect('equal')
    if i == 2:
        ax03.set_xlabel('Solar X (arcsec)')
    ax03.set_ylabel('Solar Y (arcsec)')
    ax03.set_title('IRIS SJI 2832 $\AA$')
    ax03.set_xticks(np.arange(0, 80, 20))
    ax03.set_yticks(np.arange(0, 80, 20))
    data_2832 = fits.getdata(file_2832[num])
    hdr_2832 = fits.getheader(file_2832[num])
    iris_shape = data_2832.shape
    iris_extent = hdr_2832['cdelt1']*np.array([0, iris_shape[1], 0, iris_shape[2]])
    im03 = ax03.imshow(data_2832[0], origin='lower', cmap='Purples_r', 
                       extent=iris_extent)
    im03.set_clim(-50, 550)
    t03 = ax03.text(7, 55, '  '+Time(hdr_2832['date_obs']).iso[0:-4]+' UT\n~'+
                           Time(hdr_2832['date_end']).iso[0:-4]+' UT' )
    t03.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))
    
    clump_corr_no = fits.getdata(file_corr[num]).astype(int)
    clump_max = np.max(clump_corr_no)
    oc_x = np.zeros(clump_max-1)
    oc_y = np.zeros(clump_max-1)
    for j, no in enumerate(range(1, clump_max)):
        zp, _, _ = np.where(clump_corr_no == no)
        yp, xp = np.where(clump_corr_no[zp[0]] == no)
        oc_x[j] = np.mean(xp)*hdr_2832['cdelt1']
        oc_y[j] = np.mean(yp)*hdr_2832['cdelt1']
    p03 = ax03.plot(oc_x, oc_y, ' .', alpha=0.5, color='yellow')
    t_start = Time(hdr_2832['date_obs']).jd
    t_end = Time(hdr_2832['date_end']).jd
    p020 = ax02.plot(hmi_extent[0:2], t_start*np.array([1, 1]), 
                     color='purple', alpha=0.5)
    p021 = ax02.plot(hmi_extent[0:2], t_end*np.array([1, 1]), 
                     color='purple', alpha=0.5)
    
    ax03_pos = ax03.get_position().corners()
    [ax03_xp1, ax03_yp1] = ax03_pos[0]
    [ax03_xp2, ax03_yp2] = ax03_pos[1]
    [ax02_xp1, ax02_yp1] = ax02.transData.transform((ax02.get_xlim()[1], t_end))
    [ax02_xp2, ax02_yp2] = ax02.transData.transform((ax02.get_xlim()[1], t_start))
    ax02_xp1 = ax02_xp1/width
    ax02_xp2 = ax02_xp2/width
    ax02_yp1 = ax02_yp1/height
    ax02_yp2 = ax02_yp2/height

    p001, = ax00.plot([ax02_xp1, ax03_xp1], [ax02_yp1, ax03_yp1], 
                     color='purple', alpha=0.5, ls=':')
    p002, = ax00.plot([ax02_xp2, ax03_xp2], [ax02_yp2, ax03_yp2], 
                     color='purple', alpha=0.5, ls=':')
ax00.set_axis_off()

