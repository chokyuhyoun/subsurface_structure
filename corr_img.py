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

path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)
file_corr = glob(path+r'/clump_corr*.fits')
file_2796 = glob(path+r'/fp_2796*.fits')
file_2832 = glob(path+r'/fp_2832*.fits')
num = 0
save_path = path+r'/corr_img/'+str(num)
try:
    os.mkdir(save_path)
except:
    pass

data_2796 = fits.getdata(file_2796[num])
data_2832 = fits.getdata(file_2832[num])
hdr_2796 = fits.getheader(file_2796[num])
clump_corr_no = fits.getdata(file_corr[num]).astype(int)
shape = data_2796.shape
cdelt = hdr_2796['CDELT1']
dt_2796 = TimeDelta(hdr_2796['cdelt3'], format='sec')
t_2796 = Time(hdr_2796['date_obs'])+dt_2796*np.arange(shape[0])
extent = [0, shape[2]*cdelt, 0, shape[1]*cdelt]

fig01, ax01 = plt.subplots(figsize=[6, 6])
ax01.set_xlabel('Solar X (arcsec)')
ax01.set_ylabel('Solar Y (arcsec)')
col_list = np.array(['red', 'magenta', 'yellow', 'green', 'cyan', 'blue', 'purple'])
col_list = plt.get_cmap('Paired').colors
n_col_list = len(col_list)
for i in range(shape[0]):
    print(i)
    if i == 0:
        im01 = ax01.imshow(data_2796[i], origin='lower', extent=extent, 
                           cmap='gray', alpha=0.6)
        im01.set_clim(-5, 5)
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
    for lv in no_list:
        temp_map = clump_corr_no[i] + 0
        temp_map[temp_map != lv] = 0
        dum = ax01.contour(temp_map, [lv-0.1], origin='lower', 
                           extent=extent, colors=[col_list[lv % n_col_list]], 
                           linewidths=2)
        cont_list.append(dum)

    # if len(no_list) != 0:
    #     cont02 = ax01.contour(clump_corr_no[i], no_list-0.1, origin='lower', 
    #                           extent=extent, colors=col_list[no_list % n_col_list], 
    #                           linewidths=2)
    fig01.savefig(save_path+'\\'+str(num)+'_'+f'{i:03}'+'.png', dpi=300)