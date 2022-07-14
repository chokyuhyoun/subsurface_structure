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
from scipy import interpolate
import pickle

plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
path = r'C:/Users/chokh/.spyder-py3/20151216'
os.chdir(path)
file_corr = glob(path+r'/clump_corr*.fits')
file_2796 = glob(path+r'/fp_2796*.fits')
file_2832 = glob(path+r'/fp_2832*.fits')
file_hmi = glob(path+r'/hmi_aligned*.fits')
file_vector = glob(path+r'/hmi_vector.fits')
hmi_data = fits.getdata(file_hmi[0])
hmi_vector = fits.getdata(file_vector[0])
hmi_shape = np.array(hmi_data.shape)
hdr_2832 = fits.getheader(file_2832[0])
with open(path+r'/hmi_info.p', 'rb') as file:
    hmi_time_str = pickle.load(file)
    hmi_ycxc = pickle.load(file)
    hmi_match = pickle.load(file)
with open(path+r'/hmi_vector_info.p', 'rb') as file:
    hmi_vector_time_str = pickle.load(file)
hmi_time = Time(hmi_time_str)
hmi_vector_time = Time(hmi_vector_time_str)

hmi_xp = (np.arange(hmi_shape[1])-0.5*(hmi_shape[1]-1))*0.5 + 1.64
hmi_yp = (np.arange(hmi_shape[2])-0.5*(hmi_shape[2]-1))*0.5 + 1.03
hmi_pix2data_x = interpolate.interp1d(np.arange(hmi_shape[1]), hmi_xp)
hmi_pix2data_y = interpolate.interp1d(np.arange(hmi_shape[2]), hmi_yp)
hmi_data2pix_x = interpolate.interp1d(hmi_xp, np.arange(hmi_shape[1]))
hmi_data2pix_y = interpolate.interp1d(hmi_yp, np.arange(hmi_shape[2]))

iris_xp = (np.arange(hdr_2832['naxis1'])-0.5*(hdr_2832['naxis1']-1))*hdr_2832['cdelt1']
iris_yp = (np.arange(hdr_2832['naxis2'])-0.5*(hdr_2832['naxis2']-1))*hdr_2832['cdelt2']
iris_pix2data_x = interpolate.interp1d(np.arange(hdr_2832['naxis1']), iris_xp)
iris_pix2data_y = interpolate.interp1d(np.arange(hdr_2832['naxis2']), iris_yp)

fig01, ax0n = plt.subplots(3, 2, figsize=[9, 9])
ax00 = fig01.add_axes([0, 0, 1, 1])
ax00.set_xlim(0, 1)
ax00.set_ylim(0, 1)
ax00.set_zorder(-1)
height = fig01.get_figheight()*fig01.get_dpi()
width = fig01.get_figwidth()*fig01.get_dpi()

ax01 = fig01.add_axes([0.07, 0.7, 0.25, 0.25])
ax01.set_xlabel('Solar X (arcsec)')
ax01.set_ylabel('Solar Y (arcsec)')
ax01.set_xticks(np.arange(-30, 40, 15))
ax01.set_yticks(ax01.get_xticks())
hmi_extent =  np.concatenate([hmi_xp[[0, -1]], hmi_yp[[0, -1]]])
iris_extent = np.concatenate([iris_xp[[0, -1]], iris_yp[[0, -1]]])
im01 = ax01.imshow(hmi_data[0], origin='lower', extent=hmi_extent, cmap='gray')
im01.set_clim(0, 6e4)
ax01.set_xlim(misc.minmax(iris_xp))
ax01.set_ylim(misc.minmax(iris_yp))
slit_ypos = 92
p01, = ax01.plot(ax01.get_xlim(), hmi_pix2data_y(slit_ypos)*[1, 1], 'r')
ax01.set_title('HMI Intensity')
t01 = ax01.text(ax01.get_xlim()[0]+5, ax01.get_ylim()[1]-8, hmi_time[0].iso[0:-4]+' UT')
t01.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))

ax02 = fig01.add_axes([ax01.get_position().x0, 0.07, 
                       ax01.get_position().width, 0.55])
ax02.set_xlabel('Solar X (arcsec)')
ax02.set_ylabel('Time (UT)')
ax02.set_xticks(ax01.get_xticks())
ax02.set_title('Time-Distance Map')
td_extent = np.array([hmi_extent[0:2], hmi_time.jd[[0, -1]]]).flatten()
im02 = ax02.imshow(hmi_data[:, slit_ypos, :], origin='lower', extent=td_extent, 
                   cmap='gray')
im02.set_clim(im01.get_clim())
ax02.set_aspect('auto')
td_ytick_jd = np.fix(hmi_time[0].jd)+np.array([1.5, 2.5, 3.5])
td_yticklabel = [iso[0:10]+'\n00:00:00' \
                 for iso in Time(td_ytick_jd, format='jd').iso]
ax02.set_yticks(td_ytick_jd)
ax02.set_yticklabels(td_yticklabel, rotation=90, ma='center', va='center')
ax02.set_xlim(ax01.get_xlim())
oc_b = []
umb_b = []
oc_int = []
umb_int = []
for i, ax0n1 in enumerate(ax0n):
    ax03 = ax0n1[0]
    ax04 = ax0n1[1]
    num = i*2
    ax03.set_position([0.4, ax02.get_position().y0+i*0.315, 
                       ax01.get_position().width, ax01.get_position().width])
    ax03.set_aspect('equal')
    if i == 0:
        ax03.set_xlabel('Solar X (arcsec)')
    ax03.set_ylabel('Solar Y (arcsec)')
    ax03.set_title('IRIS SJI 2832 $\AA$')
    ax03.set_xticks(ax01.get_xticks())
    ax03.set_yticks(ax01.get_xticks())
    ax04.set_position([ax03.get_position().x1, ax02.get_position().y0+i*0.315, 
                       ax01.get_position().width, ax01.get_position().width])
    ax04.set_aspect('equal')
    if i == 0:
        ax04.set_xlabel('Solar X (arcsec)')
    ax04.set_title(r'HMI $\mathbf{B}$ Strength')
    ax04.set_xticks(ax01.get_xticks())
    ax04.set_yticks(ax01.get_xticks())
    ax04.set_yticklabels([])
    ax04.set_xlim(misc.minmax(iris_xp))
    ax04.set_ylim(misc.minmax(iris_yp))
    
    data_2832 = fits.getdata(file_2832[num])
    hdr_2832 = fits.getheader(file_2832[num])
    iris_shape = data_2832.shape
    cdelt = hdr_2832['CDELT1']
    iris_dt = TimeDelta(hdr_2832['cdelt3'], format='sec')
    iris_time = Time(hdr_2832['date_obs'])+iris_dt*np.arange(iris_shape[0])
    im03 = ax03.imshow(data_2832[0], origin='lower', cmap='Purples_r', 
                       extent=iris_extent)
    im03.set_clim(-50, 700)
    t03 = ax03.text(ax03.get_xlim()[0]+5, ax03.get_ylim()[1]-12, 
                    '  '+Time(hdr_2832['date_obs']).iso[0:-4]+' UT\n~'+
                         Time(hdr_2832['date_end']).iso[0:-4]+' UT' )
    t03.set_bbox(dict(facecolor='white', alpha=0.5))
   
    with open(path+r'\flag_set_'+str(num)+'.p', 'rb') as file:
        clump_no = pickle.load(file)
        flag_set = pickle.load(file)
        flag_end = pickle.load(file)
    clump_init = np.where(flag_end == 2)[0]
    total_osc_cen = len(clump_init)
    oc_iris_pix_x = np.zeros(total_osc_cen)
    oc_iris_pix_y = np.zeros(total_osc_cen)
    oc_iris_pix_z = np.zeros(total_osc_cen)
    oc_int_part = []
    oc_b_part = []
    oc_data_x = []
    oc_data_y = []
    for j, no in enumerate(clump_init):
        zp, yp, xp = np.where(clump_no == no)  #222222222222222222222222
        oc_iris_pix_x[j] = np.mean(xp)
        oc_iris_pix_y[j] = np.mean(yp)
        oc_int_interp = interpolate.interp2d(np.arange(0, 200), np.arange(0, 200), 
                                              data_2832[zp[0]])
        oc_int_part.append(oc_int_interp(oc_iris_pix_y[j], oc_iris_pix_x[j])[0])
        
        oc_data_x.append(iris_pix2data_x(oc_iris_pix_x[j]))
        oc_data_y.append(iris_pix2data_y(oc_iris_pix_y[j]))
        vector_match = misc.arr_eq(hmi_vector_time.jd, iris_time[zp[0]].jd) #@@@@@@@@@@@@@@@@@@@@@
        oc_b_interp = interpolate.interp2d(hmi_yp, hmi_xp, hmi_vector[vector_match])    
        oc_b_part.append(oc_b_interp(oc_data_y[j], oc_data_x[j])[0])
       
    oc_int.append(oc_int_part)    
    oc_b.append(oc_b_part)
    
    p03 = ax03.plot(oc_data_x, oc_data_y, 'o', ms=3, mew=0,
                    alpha=0.5, color='yellow')
    t_start = Time(hdr_2832['date_obs']).jd
    t_end = Time(hdr_2832['date_end']).jd
    t031 = ax03.text(ax03.get_xlim()[1]-3, ax03.get_ylim()[0]+3, 
                     'N = '+str(total_osc_cen), ha='right', va='bottom')

    if i == 0 :
        leg03 = ax03.legend(p03, ['Oscillation \nCenter'], loc='lower left')
        leg03.get_frame().set_color('silver')
        leg03.get_frame().set_alpha(0.5)
    p020 = ax02.plot(hmi_extent[0:2], t_start*np.array([1, 1]), 
                     color='slateblue', alpha=0.5)
    p021 = ax02.plot(hmi_extent[0:2], t_end*np.array([1, 1]), 
                     color='slateblue', alpha=0.5)
    
    ax03_pos = ax03.get_position().corners()
    [ax03_xp1, ax03_yp1] = ax03_pos[0]
    [ax03_xp2, ax03_yp2] = ax03_pos[1]
    [ax02_xp1, ax02_yp1] = ax02.transData.transform((ax02.get_xlim()[1], t_start))
    [ax02_xp2, ax02_yp2] = ax02.transData.transform((ax02.get_xlim()[1], t_end))
    ax02_xp1 = ax02_xp1/width
    ax02_xp2 = ax02_xp2/width
    ax02_yp1 = ax02_yp1/height
    ax02_yp2 = ax02_yp2/height

    p001, = ax00.plot([ax02_xp1, ax03_xp1], [ax02_yp1, ax03_yp1], 
                     color='slateblue', alpha=0.5, ls=':')
    p002, = ax00.plot([ax02_xp2, ax03_xp2], [ax02_yp2, ax03_yp2], 
                     color='slateblue', alpha=0.5, ls=':')
    
    hmi_match2 = misc.arr_eq(hmi_vector_time.jd, iris_time[0].jd)
    im04 = ax04.imshow(hmi_vector[hmi_match2], origin='lower', extent=hmi_extent, 
                       cmap='gnuplot')
    im04.set_clim(0, 3e3)
    cont04 = ax04.contour(data_2832[0], [50, 300], origin='lower', 
                          extent=iris_extent, colors='k', linewidths=0.5)
    t04 = ax04.text(ax04.get_xlim()[0]+5, ax04.get_ylim()[1]-8, 
                    hmi_vector_time[i].iso[0:-4]+' UT', color='white')

    cax04 = fig01.add_axes([ax04.get_position().x1, ax04.get_position().y0, 
                            0.01, ax04.get_position().height])
    cb04 = fig01.colorbar(im04, cax = cax04)
    cb04.ax.set_ylabel(r'$\mathbf{B}$ strength (G)')

    sz = data_2832.shape
    cenx, ceny = 100, 95   # image coord. != array
    r1, r2 = 45, 35     # major(x) minor(y) axis
    zzp, yyp, xxp = np.mgrid[0:sz[0], 0:sz[1], 0:sz[2]]
    ell_surf_iris = ((xxp-cenx)/r1)**2+((yyp-ceny)/r2)**2 
    umb_fil01 = ell_surf_iris <= 1
    umb_int.append(data_2832[umb_fil01])
    
    cenx_data = iris_pix2data_x(cenx)
    ceny_data = iris_pix2data_y(ceny)
    r1_data = r1*hdr_2832['CDELT1']
    r2_data = r2*hdr_2832['CDELT2']
    xxp_hmi_data, yyp_hmi_data = np.meshgrid(hmi_xp, hmi_yp)
    ell_surf_hmi = ((xxp_hmi_data-cenx_data)/r1_data)**2+((yyp_hmi_data-ceny_data)/r2_data)**2 
    umb_fil02 = ell_surf_hmi <= 1
    umb_b.append(hmi_vector[hmi_match2][umb_fil02])
    
    p031 = ax03.contour(iris_xp, iris_yp, ell_surf_iris[0], [1.], 
                         origin='lower', colors='white', alpha=0.4)
    p041 = ax04.contour(hmi_xp, hmi_yp, ell_surf_hmi, [1.], 
                         origin='lower', colors='white', alpha=0.4)
    
ax00.set_axis_off()
fig01.savefig('oscillation centers.pdf', dpi=300)
#%%

fig2, ax2 = plt.subplots(3, 1, figsize=[5, 7], sharex=True, gridspec_kw=dict(hspace=0))
count2, bin2, hist2 = [], [], []
count2_umb, bin2_umb, hist2_umb = [], [], []
bins_int = np.arange(5, 50, 1)
for j, ax2n in enumerate(ax2):
    i = 2 - j
    ax2n.set_xlim(5, 30)
    ax2n.set_ylim(0, 18)
    count20, bin20, hist20 = ax2n.hist(umb_int[i], bins=bins_int+0.2, 
                                       weights=np.ones_like(umb_int[i])*100./len(umb_int[i]), 
                                       rwidth=0.5, color='lightgray', ec='black')
    count21, bin21, hist21 = ax2n.hist(oc_int[i], bins=bins_int, 
                                       weights=np.ones_like(oc_int[i])*100./len(oc_int[i]), 
                                       rwidth=0.5, color='black', ec='black')
    count2.append(count20)
    count2_umb.append(count21)
ax2[2].set_xlabel(r'IRIS SJI 2832 $\AA$ Intensity (DN)')
ax2[2].set_ylabel('Proportion (%)')
fig2.tight_layout()  

fig3, ax3 = plt.subplots(3, 1, figsize=[5, 7], sharex=True, gridspec_kw=dict(hspace=0))
count3, bin3, hist3 = [], [], []
count3_umb, bin3_umb, hist3_umb = [], [], []
bins_b = np.arange(1400, 3000, 50)
for j, ax3n in enumerate(ax3):
    i = 2 - j
    ax3n.set_xlim(1400, 3000)
    ax3n.set_ylim(0, 18)
    count30, bin30, hist30 = ax3n.hist(umb_b[i], bins=bins_b+10, 
                                       weights=np.ones_like(umb_b[i])*100./len(umb_b[i]), 
                                       rwidth=0.5, color='lightgray', ec='black')
    count31, bin31, hist31 = ax3n.hist(oc_b[i], bins=bins_b, 
                                       weights=np.ones_like(oc_b[i])*100./len(oc_b[i]), 
                                       rwidth=0.5, color='black', ec='black')
    count3.append(count30)
    count3_umb.append(count31)
ax3[2].set_xlabel(r'HMI $\mathbf{B}$ strength (G)')
ax3[2].set_ylabel('Proportion (%)')
fig3.tight_layout()  

# ax21.set_xlabel(r'IRIS SJI 2832 $\AA$ Intensity')
# ax21.set_ylabel(r'# of Oscillation Centers')
# ax21.set_xlim(5, 50)
# col = ['r', 'g', 'b']
# count20, bin20, hist20 = ax20.hist(oc_b, bins=np.arange(1400, 3000, 100), 
#                                    histtype='bar', color=col, 
#                                    alpha=1, fill=True)
# count21, bin21, hist21 = ax21.hist(oc_int, bins=np.arange(5, 300, 5), histtype='bar', color=col, 
#                                    alpha=1, fill=True)

# for i in range(3):
#     t2 = ax20.text(ax20.get_xlim()[0]+50, ax20.get_ylim()[1]-4-4*i, 
#                    hmi_vector_time[i].iso[0:-4]+' UT', color=col[i])
# fig2.tight_layout()  

# fig2.savefig('B strength hist.pdf', dpi=300)
