# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:19:04 2020

@author: chokh
"""
import os
import numpy as np
# import astropy.units as u
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.colors
from astropy.io import fits
from astropy.time import Time, TimeDelta
from mylib import nave, image_process, clump_find, misc
import scipy.ndimage as spnd
import pickle
from scipy import interpolate

path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)
num = 0

file_2796 = glob(path+r'\fp_2796*.fits')
file_2832 = glob(path+r'\fp_2832*.fits')
data_2796 = fits.getdata(file_2796[num])
data_2832 = fits.getdata(file_2832[num])
hdr_2796 = fits.getheader(file_2796[num])
sz = np.array(data_2796.shape)

cenx, ceny = 100, 95   # image coord. != array
r1, r2 = 45, 35     # major(x) minor(y) axis
yyp, xxp = np.mgrid[0:sz[1], 0:sz[2]]
dum01 = ((xxp-cenx)/r1)**2+((yyp-ceny)/r2)**2 
umb_fil = dum01 <= 1

cal = 0
  
clump_no = np.zeros(sz, dtype=int)        
c_num = 0
# f_image = np.zeros(sz)

if cal:
    for i in range(sz[0]):
        image = spnd.gaussian_filter(data_2796[i, :, :]*umb_fil, 1)
        res = clump_find.clump_find(image, 0, np.pi*3**2, dist_sq_crit=1) 
        # Wittmann (1969) D = 2700 km -> radius = 2 arcsec -> 6 pix
        # if i == 235: import pdb; pdb.set_trace()
        nothing = res == 0
        res += c_num
        res[nothing] = 0
        clump_no[i, :, :] = res
        c_num = np.max(res)
    
    last_clump = np.max(clump_no[-2])
    base_clump_num = np.max(clump_no[0, :, :])
    link_max_period = 20 
    link_max_simult = 5
    link = np.zeros([link_max_period, link_max_simult])[None, :, :] - 1 # [link, time, simult]

    for i in range(base_clump_num+1, last_clump+1):
        if i not in link:
            link = np.append(link, [np.zeros([link_max_period, link_max_simult])-1], axis=0)
            link[-1, 0, 0] = i
        zp, yp, xp = np.where(clump_no == i)
        before_temp = clump_no[zp[0]]*0
        before_temp[clump_no[zp[0]] == i] = 1
        after_temp = clump_no[zp[0]+1]
        minus = after_temp*2 - before_temp
        value, counts = np.unique(minus, return_counts=True)
        match = (value % 2 == 1) & (value > 0)
        
        if np.sum(match) != 0:
            j = ((value[match]+1)*0.5).astype(int)
            pos = np.where(link == i)
            empty_pos = np.min(np.where(link[pos[0], pos[1]+1] == -1))
            link[pos[0], pos[1]+1, empty_pos:empty_pos+len(j)] = j
 
    crit_t = 5
    mem_no = np.sum(link[:, :, 0] != -1, axis=1)
    real_link = link[mem_no >= crit_t]
    
    with open(path+r'\\link_'+str(num)+'.p', 'wb') as file:
        pickle.dump(clump_no, file)
        pickle.dump(link, file)
        pickle.dump(real_link, file)
else:
    with open(path+r'\\link_'+str(num)+'.p', 'rb') as file:
        clump_no = pickle.load(file)
        link = pickle.load(file)
        real_link = pickle.load(file)

#%%  make clump correlation movie
if 1:

    c1 = np.sqrt(np.arange(256))*np.sqrt(255.)
    c2 = np.arange(256)**2/255.    
    c3 = (c1+c2/2.)*255./(np.max(c1) + np.max(c2)/2.)
    cmap_2796 = matplotlib.colors.ListedColormap(np.array([c1, c3, c2]).T/255.)
    
    save_path = path+r'/corr_img02'
    try:
        os.mkdir(save_path)
    except:
        pass
    shape = data_2796.shape
    cdelt = hdr_2796['CDELT1']
    dt_2796 = TimeDelta(hdr_2796['cdelt3'], format='sec')
    t_2796 = Time(hdr_2796['date_obs'])+dt_2796*np.arange(shape[0])
    iris_xp = (np.arange(hdr_2796['naxis1'])-0.5*(hdr_2796['naxis1']-1))*hdr_2796['cdelt1']
    iris_yp = (np.arange(hdr_2796['naxis2'])-0.5*(hdr_2796['naxis2']-1))*hdr_2796['cdelt2']
    iris_pix2data_x = interpolate.interp1d(np.arange(hdr_2796['naxis1']), iris_xp)
    iris_pix2data_y = interpolate.interp1d(np.arange(hdr_2796['naxis2']), iris_yp)
    extent = np.concatenate([iris_xp[[0, -1]], iris_yp[[0, -1]]])
    
    fig01, ax01 = plt.subplots(figsize=[6, 6])
    ax01.set_xlabel('Solar X (arcsec)')
    ax01.set_ylabel('Solar Y (arcsec)')
    ax01.set_xlim(-20, 20)
    ax01.set_ylim(-20, 20)
    col_list = np.array(['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:pink', 'tab:cyan', \
                         'teal'])
    # col_list = np.array(plt.get_cmap('Paired').colors)
    n_col_list = len(col_list)
    for i in range(shape[0]):
        if i == 0:
            im01 = ax01.imshow(data_2796[i], origin='lower', extent=extent, 
                               alpha=0.7, cmap=cmap_2796)
            im01.set_clim(-5, 5)
        else:        
            im01.set_data(data_2796[i])
        try:
            for dum in cont01.collections:
                dum.remove()
            for dum in cont_list:
                for dum1 in dum.collections:
                    dum1.remove()
            for dum in cenplot_list:
                dum.remove()
        except:
            pass
        ax01.set_title('IRIS SJI 2796 $\AA$ '+t_2796[i].iso[0:-4]+' UT')
        cont01 = ax01.contour(data_2832[i], [100, 350], origin='lower', 
                              extent=extent, linewidths=1, colors='black')
        no_list = np.unique(clump_no[i])[1:]
        cont_list = []
        cenplot_list = []
        for lv in no_list:
            temp_map = clump_no[i] + 0 
            temp_map[temp_map != lv] = 0
            event_no = np.where(real_link == lv)[0]
            print(lv, event_no)
            if len(event_no) == 0:
                dum = ax01.contour(temp_map, levels=[lv-0.1], origin='lower', 
                                   extent=extent, colors='gray', 
                                   linestyles='dotted')
                cont_list.append(dum)
            else:
                col = col_list[event_no % n_col_list]
                dum = ax01.contour(temp_map, levels=np.linspace(0, lv-1, len(event_no)+1)[1:], 
                                   origin='lower', 
                                   extent=extent, colors=col)

                cont_list.append(dum)
                if len(event_no) == 1:
                    if real_link[event_no[0], 0, 0] == lv:
                        lv_pos = np.where(clump_no == lv)
                        xcen = np.mean(lv_pos[2])
                        ycen = np.mean(lv_pos[1])
                        xcen_data = iris_pix2data_x(xcen)
                        ycen_data = iris_pix2data_y(ycen)
                        p02, = ax01.plot(xcen_data, ycen_data, '+', 
                                         ms=10, mew=2, color=col[0])
                        cenplot_list.append(p02)
            # import pdb; pdb.set_trace()
        fig01.savefig(save_path + r'/' + f'{i:03}'+'.png', dpi=300)
        
    img_list = glob(save_path+'*.png')
    res = misc.ffmpeg(img_list, 20)
    
 



#%%



        
flag_num = 1
crit_t = 5
clump_corr_no = np.zeros(clump_no.shape)
clump_corr_init = np.zeros(clump_no.shape)
for i in set(flag_set):
    if i == 0:
        continue
    if i in clump_no[0]:
        continue
    clump_number = np.where(flag_set == i)[0]
    start_t = np.where(clump_no == clump_number[0])[0][0]
    end_t = np.where(clump_no == clump_number[-1])[0][0]
    if (end_t - start_t >= crit_t) & (start_t != 0):
        # print(i, clump_number)
        for j in clump_number:
            clump_corr_no[clump_no == j] = flag_num
            if flag_end[j] == 1:
                # import pdb; pdb.set_trace()
                flag_end[j] = 2
        flag_num += 1

#%%
# clump_corr_no[clump_corr_no == 46] = 0  # manually remove
filename_2796 = 'clump_corr_'+file_2796[num][-20:]
try:
    os.remove(filename_2796)
except:
    pass    
fits.writeto(filename_2796, clump_corr_no)        

with open(path+r'\\flag_set_'+str(num)+'.p', 'wb') as file:
    pickle.dump(clump_no, file)
    pickle.dump(flag_set, file)
    pickle.dump(flag_end, file)
    
    