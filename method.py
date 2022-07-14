# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 20:21:35 2020

@author: chokh
"""
import os
import numpy as np
# import astropy.units as u
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects
import matplotlib as mpl
from matplotlib.colors import LightSource
from matplotlib.colors import colorConverter
import mpl_toolkits.mplot3d as m3d
# from mpl_toolkits.mplot3d import proj3d
from astropy.io import fits
from astropy.time import Time, TimeDelta
from mylib import misc
import scipy.ndimage as spnd
import pickle


mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

path = r'C:\Users\chokh\.spyder-py3\20151216'
os.chdir(path)
num = 0

file_2796 = glob(path+r'\fp_2796*.fits')
data_2796 = fits.getdata(file_2796[num])
hdr_2796 = fits.getheader(file_2796[num])
shape = data_2796.shape
cdelt = hdr_2796['CDELT1']
dt_2796 = TimeDelta(hdr_2796['cdelt3'], format='sec')
t_2796 = Time(hdr_2796['date_obs'])+dt_2796*np.arange(shape[0])

file_clump = glob(path+r'\clump_corr*.fits')
clump_corr_no = fits.getdata(file_clump[num])

fig01 = plt.figure(figsize=[8, 8])

# ---- observation data
ax01 = fig01.add_axes([0.1, 0.65, 0.25, 0.25])
exam_t = 152
hw = 10
fw = hw*2 + 1
exam_no = 33
zpos, ypos, xpos = np.where(clump_corr_no == exam_no)
exam_xc = (np.sum(misc.minmax(xpos))*0.5).astype(int)
exam_yc = (np.sum(misc.minmax(ypos))*0.5).astype(int)
exam_zc = (np.sum(misc.minmax(zpos))*0.5).astype(int)
part_2796 = data_2796[exam_t, exam_yc-hw:exam_yc+hw, exam_xc-hw:exam_xc+hw]
im01 = ax01.imshow(part_2796, origin='lower')
im01.set_clim(-10, 10)
ax01.set_ylabel('Solar Y (pixel)')
ax01.set_xlabel('Solar X (pixel)')
ax01.set_title('IRIS SJI 2796 $\AA$, $I(x, y)$')

# ---- Gaussian smoothing 
dum = ax01.get_position()
ax02 = fig01.add_axes([dum.x0, dum.y0-0.55, dum.width, dum.height])
part_2796_sm = spnd.gaussian_filter(part_2796, 1)
im02 = ax02.imshow(part_2796_sm, origin='lower')
im02.set_clim(im01.get_clim())
ax02.set_ylabel('Solar Y (pixel)')
ax02.set_xlabel('Solar X (pixel)')
ax02.set_title('$F(x, y)$')

x0, y0 = (hw, hw)
p021, = ax02.plot([x0], [y0], 'o', color='black', markersize=5)

ar_len = 10
grady, gradx = np.gradient(part_2796_sm)
ar021 = mpatches.FancyArrowPatch(
                    (x0-gradx[x0, y0]*ar_len*0.25, y0-grady[x0, y0]*ar_len*0.25), 
                    (x0+gradx[x0, y0]*ar_len*0.25, y0+grady[x0, y0]*ar_len*0.25))
ax02.add_patch(ar021)
ar021.set_arrowstyle('simple', head_length=5, head_width=5)
t021 = ax02.text(x0+gradx[x0, y0]*ar_len*0.2, y0+grady[x0, y0]*ar_len*0.4, 
                 r'$\mathbf{v} = \mathbf{\nabla} F$', va='center', ha='left', 
                 fontsize=15)
t021.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

grad_orthox = -grady
grad_orthoy = gradx
ar022 = mpatches.FancyArrowPatch(
                    (x0-grad_orthox[x0, y0]*ar_len*0.25, y0-grad_orthoy[x0, y0]*ar_len*0.25), 
                    (x0+grad_orthox[x0, y0]*ar_len*0.25, y0+grad_orthoy[x0, y0]*ar_len*0.25))
ax02.add_patch(ar022)
ar022.set_arrowstyle('simple', head_length=5, head_width=5)
t022 = ax02.text(x0+grad_orthox[x0, y0]*ar_len*0.28, y0+grad_orthoy[x0, y0]*ar_len*0.28, 
                 r'$\mathbf{v}_{\bot}$ ', va='bottom', ha='left', fontsize=15)
t022.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])
# t023 = fig01.text(0.11, ax02.get_position().y0-0.1, r'$\vec v = \mathbf{\nabla} F(x, y)$')

con12 = mpatches.ConnectionPatch((4, fw-4), (4, 4),
                              coordsA="data", coordsB="data",
                              arrowstyle="<|-", connectionstyle='arc3,rad=-0.3', 
                              mutation_scale=20, clip_on=False,
                              axesA=ax02, axesB=ax01)
ax02.add_artist(con12)
t12 = fig01.text(0.05, 0.5*(ax01.get_position().y0 + ax02.get_position().y1), 
                 r'${F(x, y)} = I(x, y)*G(\sigma)$', va='center', fontsize=15, 
                 bbox=dict(boxstyle='round', facecolor='white'))

# ---- Binary map
dum = ax02.get_position()
ax03 = fig01.add_axes([dum.x0+0.6, dum.y0, dum.width, dum.height], zorder=10)
part_clump = clump_corr_no[exam_t, exam_yc-hw:exam_yc+hw, exam_xc-hw:exam_xc+hw]
cmap03 = mpl.colors.LinearSegmentedColormap.from_list("", ["white","red"])
part_clump1 = part_clump + 0
part_clump1[part_clump != (exam_no-1)] = np.nan
im031 = ax03.imshow(part_clump1, origin='lower', cmap='binary', alpha=0.2)
part_clump2 = part_clump + 0
part_clump2[part_clump != (exam_no)] = np.nan
im032 = ax03.imshow(part_clump2, origin='lower', cmap=cmap03)
im031.set_clim(0, 2)
im032.set_clim(0, 2)
ax03.set_ylabel('Solar Y (pixel)')
ax03.set_xlabel('Solar X (pixel)')
ax03.set_title('$B(x, y)$')

con23 = mpatches.ConnectionPatch((4, 4), (fw-4, 4),
                              coordsA="data", coordsB="data",
                              arrowstyle="<|-", connectionstyle='arc3,rad=-0.3', 
                              mutation_scale=20, clip_on=False,
                              axesA=ax03, axesB=ax02)
ax03.add_artist(con23)
t23_str = r'$B(x, y)  $'+'\n' \
        + r'  $= 1$  if  $\frac{\partial^2 F}{\partial \mathbf{v}^2} \leq 0 $'+'\n' \
        + r'          $&\/\/ \frac{\partial^2 F}{\partial \mathbf{v}_{\bot}^2} \leq 0$'+'\n' \
        + r'  $= 0$   otherwise'
t23 = fig01.text(0.4, 0.15,
                 t23_str, va='center', ha='left', fontsize=15, 
                 bbox=dict(boxstyle='round', facecolor='white'), zorder=10)
t03 = ax03.text(1, 18, r'Area '+'\n'+r' $\geq$ $3^2\pi$ pixel'+'\n'+r' $\sim \pi$ arcsec$^2$', 
                color='red', va='top')


# ---- 3D patch

ax04 = fig01.add_subplot(projection='3d', position=[0.45, 0.5, 0.45, 0.45])
ax04.set_xticks(np.arange(0, fw, 5))
ax04.set_yticks(np.arange(0, fw, 5))
ax04.set_zticks(np.arange(0, fw, 5))
ax04.set_xlim(ax01.get_xlim())
ax04.set_ylim(ax01.get_ylim())
ax04.set_zlim(-0.5, np.max(zpos)-np.min(zpos)+1.5)
ax04.set_xlabel('Solar X (pixel)')
ax04.set_ylabel('Solar Y (pixel)')
ax04.set_zlabel('Time (pixel)')
ax04.view_init(34, -73)

axmin = ax01.get_xlim()[0]
axmax = ax01.get_xlim()[1]
verts = np.array([[[axmin, axmin], [axmin, axmax], [axmax, axmax], \
                   [axmax, axmin], [axmin, axmin]]])
im04_zlev = exam_t-np.min(zpos)
im04 = mpl.collections.PolyCollection(verts, facecolors='gray', alpha=0.6)
ax04.add_collection3d(im04, zs=im04_zlev, zdir='z')
# yyp, xxp = misc.grid(part_2796)
# zzp = xxp*0 + im04_zlev
# ax04.plot_surface(xxp, yyp, zzp, color=colorConverter.to_rgba('gray', alpha=0.6), zorder=-1)

part_vol = clump_corr_no[np.min(zpos):np.max(zpos)+1, exam_yc-hw:exam_yc+hw, exam_xc-hw:exam_xc+hw]
part_vol[part_vol != exam_no] = 0
fc = np.zeros(part_vol.shape).astype(str)
# fc[part_vol != 0] = '#00000005'
fc[part_vol == exam_no] = '#FF0000FF'
zp, yp, xp = np.indices(np.array(part_vol.shape)+1)
v04 = ax04.voxels(xp, yp, zp, (part_vol != 0), facecolors=fc, edgecolors='#80808030', 
                  linewidth=0.2, shade=True, lightsource=LightSource(azdeg=270, altdeg=-45))

ar04 = misc.Arrow3D([3, 3], [-1, -1], [0, ax04.get_zlim()[1]], mutation_scale=20, 
                    arrowstyle='<|-|>', color='k', zorder=10)
ax04.add_artist(ar04)
t04 = ax04.text(-6, 0, im04_zlev+1, 
                r'Duration '+'\n'+r'  $\geq$ 5 pixel'+'\n'+r'  $\sim$ 60 sec', 
                color='k', va='bottom', ha='left', zorder=10)

# ---- connecting lines
[ax03_xp1, ax03_yp1] = ax03.transData.transform((ax03.get_xlim()[0], ax03.get_ylim()[0]))
[ax03_xp2, ax03_yp2] = ax03.transData.transform((ax03.get_xlim()[1], ax03.get_ylim()[1]))
dumx, dumy, _ = m3d.proj3d.proj_transform(ax03.get_xlim(), ax03.get_ylim(), 
                                      im04_zlev*np.array([1, 1]), ax04.get_proj())
[ax04_xp1, ax04_yp1] = ax04.transData.transform((dumx[0], dumy[0]))
[ax04_xp2, ax04_yp2] = ax04.transData.transform((dumx[1], dumy[1]))

ax00 = fig01.add_axes([0, 0, 1, 1])
ax00.set_xlim(0, 1)
ax00.set_ylim(0, 1)
height = fig01.get_figheight()*fig01.get_dpi()
width = fig01.get_figwidth()*fig01.get_dpi()
ax03_xp1 = ax03_xp1/width
ax03_xp2 = ax03_xp2/width
ax04_xp1 = ax04_xp1/width
ax04_xp2 = ax04_xp2/width
ax03_yp1 = ax03_yp1/height
ax03_yp2 = ax03_yp2/height
ax04_yp1 = ax04_yp1/height
ax04_yp2 = ax04_yp2/height

transFigure = fig01.transFigure.inverted()
coord031 = transFigure.transform(ax00.transData.transform([ax03_xp1, ax03_yp1]))
coord032 = transFigure.transform(ax00.transData.transform([ax03_xp2, ax03_yp2]))
coord041 = transFigure.transform(ax00.transData.transform([ax04_xp1, ax04_yp1]))
coord042 = transFigure.transform(ax00.transData.transform([ax04_xp2, ax04_yp2]))

line01 = mpl.lines.Line2D((coord031[0], coord041[0]), (coord031[1], coord041[1]), 
                          transform=fig01.transFigure, color='gray')
line02 = mpl.lines.Line2D((coord032[0], coord042[0]), (coord032[1], coord042[1]), 
                          transform=fig01.transFigure, color='gray')
fig01.lines = [line01, line02]
ax00.set_visible(False)

fig01.set_rasterized(True)
fig01.savefig('method.pdf', dpi=300)



