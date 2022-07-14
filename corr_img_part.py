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
from astropy.io import fits
from astropy.time import Time, TimeDelta
from mylib import nave, image_process, clump_find
from mylib import misc
import scipy.ndimage as spnd
from scipy import interpolate

path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)
file_corr = glob(path+r'/clump_corr*.fits')
file_2796 = glob(path+r'/fp_2796*.fits')
file_2832 = glob(path+r'/fp_2832*.fits')
num = 0
save_path = path+r'/corr_img_part/'
try:
    os.mkdir(save_path)
except:
    pass

data_2796 = fits.getdata(file_2796[num])
data_2832 = fits.getdata(file_2832[num])
hdr_2796 = fits.getheader(file_2796[num])
hdr_2832 = fits.getheader(file_2796[num])
clump_corr_no = fits.getdata(file_corr[num]).astype(int)
shape = data_2796.shape
cdelt = hdr_2796['CDELT1']
dt_2796 = TimeDelta(hdr_2796['cdelt3'], format='sec')
t_2796 = Time(hdr_2796['date_obs'])+dt_2796*np.arange(shape[0])
extent = [0, shape[2]*cdelt, 0, shape[1]*cdelt]

iris_xp = (np.arange(hdr_2832['naxis1']))*hdr_2832['cdelt1']
iris_yp = (np.arange(hdr_2832['naxis2']))*hdr_2832['cdelt2']
iris_pix2data_x = interpolate.interp1d(np.arange(hdr_2832['naxis1']), iris_xp)
iris_pix2data_y = interpolate.interp1d(np.arange(hdr_2832['naxis2']), iris_yp)


fig01, ax01 = plt.subplots(figsize=[6, 6])
ax01.set_xlabel('Solar X (arcsec)')
ax01.set_ylabel('Solar Y (arcsec)')
ax01.set_xlim(32, 53)
ax01.set_ylim(22, 43)
col_list = np.array(['red', 'magenta', 'yellow', 'green', 'cyan', 'blue', 'purple'])
col_list = plt.get_cmap('Paired').colors
n_col_list = len(col_list)

begin = 1
for i in range(168, 186):    
# for i in range(170, 171):
    print(i)
    if begin == 1:
        im01 = ax01.imshow(data_2796[i], origin='lower', extent=extent, 
                           cmap='gray', alpha=0.6)
        im01.set_clim(-5, 5)
        begin = 0
    else:        
        im01.set_data(data_2796[i])
    try:
        for dum in cont01.collections:
            dum.remove()
        for dum in cont_list:
            for dum1 in dum.collections:
                dum1.remove()
    except:
        pass
    ax01.set_title('IRIS SJI 2796 $\AA$ '+t_2796[i].iso[0:-4]+' UT')
    cont01 = ax01.contour(data_2832[i], [100, 350], origin='lower', 
                          extent=extent,linewidths=1.5, colors='black')
    no_list = np.unique(clump_corr_no[i])[1:]
    cont_list = []
    # for lv in no_list:
    lv = 25
    temp_map = clump_corr_no[i] + 0
    temp_map[temp_map != lv] = 0
    dum = ax01.contour(temp_map, [lv-0.1], origin='lower', 
                       extent=extent, colors=[col_list[lv % n_col_list]], 
                       linewidths=2)
    cont_list.append(dum)
    if i == 168:
        temp_map = clump_corr_no[170]
        temp_map[temp_map != lv] = 0
        ypos, xpos = np.where(temp_map != 0)
        ocx = iris_pix2data_x(np.mean(xpos))
        ocy = iris_pix2data_y(np.mean(ypos))
        p01, = ax01.plot(ocx, ocy, '+r', markersize=15, markeredgewidth=3)
        an01 = ax01.annotate('Oscillation Center', xy=(ocx-0.2, ocy+0.2), 
                            xytext=(38, 37.5), fontsize=13,
                            arrowprops=dict(arrowstyle='-|>', lw=2, fc='k', ec='k'))
    # if len(no_list) != 0:
    #     cont02 = ax01.contour(clump_corr_no[i], no_list-0.1, origin='lower', 
    #                           extent=extent, colors=col_list[no_list % n_col_list], 
    #                           linewidths=2)
    fig01.savefig(save_path+f'{i:03}'+'.png', dpi=300)

pngfile = glob(save_path+'???.png')
os.chdir(save_path)
dum = misc.ffmpeg(pngfile, 4)
