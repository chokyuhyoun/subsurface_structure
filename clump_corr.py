# -*- coding: utf-8 -*-
"""
Created on Thu May 28 12:19:04 2020

@author: chokh
"""
import os
import numpy as np
# import astropy.units as u
import glob
import matplotlib.pyplot as plt
from mylib.show3d import show3d
from astropy.io import fits
from mylib import nave, image_process, clump_find, misc
import scipy.ndimage as spnd
import pickle

path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)
num = 0

file_2796 = glob.glob(path+r'\fp_2796*.fits')
file_2832 = glob.glob(path+r'\fp_2832*.fits')
data_2796 = fits.getdata(file_2796[num])
data_2832 = fits.getdata(file_2832[num])
sz = np.array(data_2796.shape)

cenx, ceny = 100, 95   # image coord. != array
r1, r2 = 45, 35     # major(x) minor(y) axis
yyp, xxp = np.mgrid[0:sz[1], 0:sz[2]]
dum01 = ((xxp-cenx)/r1)**2+((yyp-ceny)/r2)**2 
umb_fil = dum01 <= 1

sav_name = os.path.basename(file_2796[num])[0:-5]+'_nave05.p'  # 25 -> 30
cal = True

   
clump_no = np.zeros(sz, dtype=int)        
c_num = 0
f_image = np.zeros(sz)
for i in range(sz[0]):
# for i in range(2):
    # print(i)
    image = spnd.gaussian_filter(data_2796[i, :, :]*umb_fil, 1)
    res = clump_find.clump_find(image, 0, np.pi*(1*3)**2, dist_sq_crit=1)
    nothing = res == 0
    res += c_num
    res[nothing] = 0
    clump_no[i, :, :] = res
    c_num = np.max(res)
    f_image[i] = image
    # print(i, mean)

total_clump = np.max(clump_no)
base_clump_num = np.max(clump_no[0, :, :])
flag = 1
crit = 0.75  # % coincidence of the area
flag_set = np.zeros(total_clump+1, dtype=int)
#%%
if cal:
    for i in range(total_clump, base_clump_num, -1):
    # for i in range(10, 11):
        # before_clump = np.zeros(sz[1:3])
        zp, yp, xp = np.where(clump_no == i)
        # area = len(xp)
        # fwhm = 2*np.sqrt(area)/np.pi
        # fwhm = np.max([np.max(xp)-np.min(xp), np.max(yp)-np.min(yp)])
        # img01 = f_image[zp[0]]
        # img02 = f_image[zp[0]-1]
        # before_x, before_y = nave.nave_track(img01, img02, fwhm, xp, yp)
        # before_clump[before_y.astype(int), before_x.astype(int)] = 1
        after_temp = clump_no[zp[0]]*0
        after_temp[clump_no[zp[0]] == i] = 1
        before_temp = clump_no[zp[0]-1]
        minus = before_temp*2 - after_temp
        value, counts = np.unique(minus, return_counts=True)
        # ratio = counts/np.sum(before_clump)
        match = (value % 2 == 1) & (value > 0)
        if np.sum(match) != 0:
            j = ((value[match]+1)*0.5).astype(int)
            if np.sum(flag_set[j]) == 0:
                if flag_set[i] == 0:
                    flag_set[j] = flag
                    flag_set[i] = flag
                    flag +=1
                else:
                    flag_set[j] = flag_set[i]
            else:
                dum = np.append(flag_set[j], flag_set[i])
                assign = min(dum[dum != 0])
                for jj in np.append(j, i):
                    if flag_set[jj] == 0:
                        flag_set[jj] = assign
                    else:
                        flag_set[flag_set == flag_set[jj]] = assign
            print(i, flag_set[i], j, flag_set[j])
        # if flag_set[0] != 0:
        # if i == 2419:
        #     import pdb; pdb.set_trace()
    with open(path+r'\\flag_set_'+str(num)+'.p', 'wb') as file:
        pickle.dump(clump_no, file)
        pickle.dump(flag_set, file)
else:
    with open(path+r'\\flag_set_'+str(num)+'.p', 'rb') as file:
        clump_no = pickle.load(file)
        flag_set = pickle.load(file)
        
flag_num = 1
crit_t = 7
clump_corr_no = np.zeros(clump_no.shape)
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
        flag_num += 1

# clump_corr_no[clump_corr_no == 46] = 0  # manually remove
filename_2796 = 'clump_corr_'+file_2796[0][-20:]
try:
    os.remove(filename_2796)
except:
    pass    
fits.writeto(filename_2796, clump_corr_no)        
    
    