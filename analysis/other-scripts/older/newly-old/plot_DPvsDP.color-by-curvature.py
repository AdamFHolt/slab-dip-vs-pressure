#!/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions_plotting import plot_BvsFullForce_wKthresh, plot_BvsDP_wKthresh
from functions_plotting import plot_BvsFullForce_Kcolored, plot_BvsDP_Kcolored

analysis_depth_dz = float(sys.argv[1])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[2])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[3])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)
depths = str(sys.argv[4])               # "normal" (230e3, 330e3, 430e3) or "other" (250e3, 300e3, 350e3)

tactual_min = 11 # first time step to use
tmin = tactual_min - 8

if depths == "normal":
    analysis_depth1 = 200e3
    analysis_depth2 = 300e3
    analysis_depth3 = 400e3
    plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.200-to-400km.K-colored.png'])
    plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.200-to-400km.K-colored.pdf'])
elif depths == "other":
    analysis_depth1 = 250e3
    analysis_depth2 = 300e3
    analysis_depth3 = 350e3
    plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.250-to-350km.K-colored.png'])
    plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.250-to-350km.K-colored.pdf'])
else:
    print("Invalid depth option. Choose 'normal' or 'other'")
    sys.exit()

name_weak_bothfree   = "2D_compositional_subd_lower-res_new_50plates"
name_weak_fixedSP    = "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name_weak_fixedOP    = "2D_compositional_subd_lower-res_new_FixedOP_50plates"
name_ref_bothfree    = "2D_compositional_subd_lower-res_new_250plates"
name_ref_fixedSP     = "2D_compositional_subd_lower-res_new_FixedSP_250plates"
name_ref_fixedOP     = "2D_compositional_subd_lower-res_new_FixedOP_250plates"
name_ref2_bothfree   = "2D_compositional_subd_lower-res_new2"         # 500
name_ref2_fixedSP    = "2D_compositional_subd_FixedSP_lower-res_new2" # 500 
name_ref2_fixedOP    = "2D_compositional_subd_FixedOP_lower-res_new" # 500
name_strong_bothfree = "2D_compositional_subd_lower-res_new_1000plates"
name_strong_fixedSP  = "2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
name_strong_fixedOP  = "2D_compositional_subd_lower-res_new_FixedOP_1000plates"
name_new_bothfree    = "2D_compositional_subd_lower-res_new_375plates"
name_new_fixedSP     = "2D_compositional_subd_lower-res_new_FixedSP_375plates"
name_new_fixedOP     = "2D_compositional_subd_lower-res_new_FixedOP_375plates"


# shall depths
text1_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# intermediate
text2_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# deep
text3_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


ref1_bothfree   = np.loadtxt((text1_ref_bothfree)) 
ref1_fixedSP    = np.loadtxt((text1_ref_fixedSP))
ref1_fixedOP    = np.loadtxt((text1_ref_fixedOP))
ref1b_bothfree   = np.loadtxt((text1_ref2_bothfree)) 
ref1b_fixedSP    = np.loadtxt((text1_ref2_fixedSP))
ref1b_fixedOP    = np.loadtxt((text1_ref2_fixedOP))
weak1_bothfree  = np.loadtxt((text1_weak_bothfree)) 
weak1_fixedSP   = np.loadtxt((text1_weak_fixedSP))
weak1_fixedOP   = np.loadtxt((text1_weak_fixedOP))
strong1_bothfree = np.loadtxt((text1_strong_bothfree)) 
strong1_fixedSP     = np.loadtxt((text1_strong_fixedSP))
strong1_fixedOP     = np.loadtxt((text1_strong_fixedOP))
new1_bothfree = np.loadtxt((text1_new_bothfree)) 
new1_fixedSP     = np.loadtxt((text1_new_fixedSP))
new1_fixedOP     = np.loadtxt((text1_new_fixedOP))

ref2_bothfree   = np.loadtxt((text2_ref_bothfree)) 
ref2_fixedSP    = np.loadtxt((text2_ref_fixedSP))
ref2_fixedOP    = np.loadtxt((text2_ref_fixedOP))
ref2b_bothfree   = np.loadtxt((text2_ref2_bothfree)) 
ref2b_fixedSP    = np.loadtxt((text2_ref2_fixedSP))
ref2b_fixedOP    = np.loadtxt((text2_ref2_fixedOP))
weak2_bothfree  = np.loadtxt((text2_weak_bothfree)) 
weak2_fixedSP   = np.loadtxt((text2_weak_fixedSP))
weak2_fixedOP   = np.loadtxt((text2_weak_fixedOP))
strong2_bothfree = np.loadtxt((text2_strong_bothfree)) 
strong2_fixedSP     = np.loadtxt((text2_strong_fixedSP))
strong2_fixedOP     = np.loadtxt((text2_strong_fixedOP))
new2_bothfree = np.loadtxt((text2_new_bothfree)) 
new2_fixedSP     = np.loadtxt((text2_new_fixedSP))
new2_fixedOP     = np.loadtxt((text2_new_fixedOP))

ref3_bothfree   = np.loadtxt((text3_ref_bothfree)) 
ref3_fixedSP    = np.loadtxt((text3_ref_fixedSP))
ref3_fixedOP    = np.loadtxt((text3_ref_fixedOP))
ref3b_bothfree   = np.loadtxt((text3_ref2_bothfree)) 
ref3b_fixedSP    = np.loadtxt((text3_ref2_fixedSP))
ref3b_fixedOP    = np.loadtxt((text3_ref2_fixedOP))
weak3_bothfree  = np.loadtxt((text3_weak_bothfree)) 
weak3_fixedSP   = np.loadtxt((text3_weak_fixedSP))
weak3_fixedOP   = np.loadtxt((text3_weak_fixedOP))
strong3_bothfree = np.loadtxt((text3_strong_bothfree)) 
strong3_fixedSP     = np.loadtxt((text3_strong_fixedSP))
strong3_fixedOP     = np.loadtxt((text3_strong_fixedOP))
new3_bothfree    = np.loadtxt((text3_new_bothfree)) 
new3_fixedSP     = np.loadtxt((text3_new_fixedSP))
new3_fixedOP     = np.loadtxt((text3_new_fixedOP))


fig=plt.figure()
gs=GridSpec(2,3) 

#### SCATTER PLOTS ####
def fixed_aspect_ratio(ratio):
    '''
    Set a fixed aspect ratio on matplotlib plots 
    regardless of axis units
    '''
    xvals,yvals = plt.gca().axes.get_xlim(),plt.gca().axes.get_ylim()

    xrange = xvals[1]-xvals[0]
    yrange = yvals[1]-yvals[0]
    plt.gca().set_aspect(ratio*(xrange/yrange), adjustable='box')


# set up color map
color_map = cm.get_cmap('inferno')
color_min = 0
color_max = 0.0025
interval  = color_max/10.
norm = matplotlib.colors.BoundaryNorm(np.arange(color_min,color_max+interval,interval), color_map.N)

# DEPTH 1
ax=fig.add_subplot(gs[0,0])
# 50
for i in range(tmin,len(weak1_bothfree)):
	im = plt.scatter(weak1_bothfree[i,4]/1.e6,(weak1_bothfree[i,3]-weak1_bothfree[i,6]+weak1_bothfree[i,17])/1.e6,s=10,c=weak1_bothfree[i,11]*1e3,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.1,zorder=2,marker='o')
plot_BvsFullForce_Kcolored(tmin,weak1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,weak1_fixedOP,'black','^',color_map,norm,zorder=2)
# 250
plot_BvsFullForce_Kcolored(tmin,ref1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,ref1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref1_fixedOP,'black','^',color_map,norm,zorder=2)
# 375
plot_BvsFullForce_Kcolored(tmin,new1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,new1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,new1_fixedOP,'black','^',color_map,norm,zorder=2)
# 500
plot_BvsFullForce_Kcolored(tmin,ref1b_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,ref1b_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref1b_fixedOP,'black','^',color_map,norm,zorder=2)
# 1000
plot_BvsFullForce_Kcolored(tmin,strong1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,strong1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,strong1_fixedOP,'black','^',color_map,norm,zorder=2)

# axis stuff
plt.xlim(-10,  37.5); plt.ylim(-10,  37.5)
plt.plot([-10, 37.5], [-10, 37.5], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel(r'$\mathregular{\Delta P +}$ slab stress  [MPa]',size=6.5)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30] )
fixed_aspect_ratio(1)

# c-bar
cbar = plt.colorbar(im, cax = fig.add_axes([0.322, 0.57, 0.008, 0.08]), orientation='vertical', ticks=[0,0.001,0.002], extend='max') # left, base, length left-right, length top-bott
cbar.ax.tick_params(axis='y',labelsize=4.5,pad=1)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.yaxis.set_label_position('left')

# DEPTH 2
ax=fig.add_subplot(gs[0,1])
# 50
plot_BvsFullForce_Kcolored(tmin,weak2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,weak2_fixedSP,'black','v',color_map,norm, zorder=3)
plot_BvsFullForce_Kcolored(tmin,weak2_fixedOP,'black','^',color_map,norm, zorder=2)
# 250
plot_BvsFullForce_Kcolored(tmin,ref2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,ref2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref2_fixedOP,'black','^',color_map,norm,zorder=2)
# 375
plot_BvsFullForce_Kcolored(tmin,new2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,new2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,new2_fixedOP,'black','^',color_map,norm,zorder=2)
# 500
plot_BvsFullForce_Kcolored(tmin,ref2b_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,ref2b_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref2b_fixedOP,'black','^',color_map,norm,zorder=2)
# 1000
plot_BvsFullForce_Kcolored(tmin,strong2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,strong2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,strong2_fixedOP,'black','^',color_map,norm,zorder=2)

# axis stuff
plt.xlim(-10,  37.5); plt.ylim(-10,  37.5)
plt.plot([-10, 37.5], [-10, 37.5], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30] )
fixed_aspect_ratio(1)

# DEPTH 3
ax=fig.add_subplot(gs[0,2])
# 50
plot_BvsFullForce_Kcolored(tmin,weak3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,weak3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,weak3_fixedOP,'black','^',color_map,norm,zorder=2)
# 250
plot_BvsFullForce_Kcolored(tmin,ref3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,ref3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref3_fixedOP,'black','^',color_map,norm,zorder=2)
# 375
plot_BvsFullForce_Kcolored(tmin,new3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,new3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,new3_fixedOP,'black','^',color_map,norm,zorder=2)
# 500
plot_BvsFullForce_Kcolored(tmin,ref3b_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,ref3b_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref3b_fixedOP,'black','^',color_map,norm,zorder=2)
# 1000
plot_BvsFullForce_Kcolored(tmin,strong3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsFullForce_Kcolored(tmin,strong3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,strong3_fixedOP,'black','^',color_map,norm,zorder=2)

# axis stuff
plt.xlim(-10,  37.5); plt.ylim(-10,  37.5)
plt.plot([-10, 37.5], [-10, 37.5], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30] )
fixed_aspect_ratio(1)

#########################################################################################
################################ Now just DP vs. DP #####################################
#########################################################################################

# DEPTH 1
ax=fig.add_subplot(gs[1,0])
# 50
plot_BvsDP_Kcolored(tmin,weak1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,weak1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,weak1_fixedOP,'black','^',color_map,norm,zorder=2)
# 250
plot_BvsDP_Kcolored(tmin,ref1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,ref1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,ref1_fixedOP,'black','^',color_map,norm,zorder=2)
# 375
plot_BvsDP_Kcolored(tmin,new1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,new1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,new1_fixedOP,'black','^',color_map,norm,zorder=2)
# 500
plot_BvsDP_Kcolored(tmin,ref1b_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,ref1b_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,ref1b_fixedOP,'black','^',color_map,norm,zorder=2)
# 1000
plot_BvsDP_Kcolored(tmin,strong1_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,strong1_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,strong1_fixedOP,'black','^',color_map,norm,zorder=2)

# axis stuff
plt.xlim(-10,  37.5); plt.ylim(-10,  37.5)
plt.plot([-10, 37.5], [-10, 37.5], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel(r'$\mathregular{\Delta P}$  [MPa]',size=6.5)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30] )
fixed_aspect_ratio(1)

# DEPTH 2
ax=fig.add_subplot(gs[1,1])
# 50
plot_BvsDP_Kcolored(tmin,weak2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,weak2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,weak2_fixedOP,'black','^',color_map,norm,zorder=2)
# 250
plot_BvsDP_Kcolored(tmin,ref2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,ref2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,ref2_fixedOP,'black','^',color_map,norm,zorder=2)
# 375
plot_BvsDP_Kcolored(tmin,new2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,new2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,new2_fixedOP,'black','^',color_map,norm,zorder=2)
# 500
plot_BvsDP_Kcolored(tmin,ref2b_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,ref2b_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,ref2b_fixedOP,'black','^',color_map,norm,zorder=2)
# 1000
plot_BvsDP_Kcolored(tmin,strong2_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,strong2_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,strong2_fixedOP,'black','^',color_map,norm,zorder=2)

# axis stuff
plt.xlim(-10,  37.5); plt.ylim(-10,  37.5)
plt.plot([-10, 37.5], [-10, 37.5], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30] )
fixed_aspect_ratio(1)

# DEPTH 3
ax=fig.add_subplot(gs[1,2])
# 50
plot_BvsDP_Kcolored(tmin,weak3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,weak3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,weak3_fixedOP,'black','^',color_map,norm,zorder=2)
# 250
plot_BvsDP_Kcolored(tmin,ref3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,ref3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,ref3_fixedOP,'black','^',color_map,norm,zorder=2)
# 375
plot_BvsDP_Kcolored(tmin,new3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,new3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,new3_fixedOP,'black','^',color_map,norm,zorder=2)
# 500
plot_BvsDP_Kcolored(tmin,ref3b_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,ref3b_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,ref3b_fixedOP,'black','^',color_map,norm,zorder=2)
# 1000
plot_BvsDP_Kcolored(tmin,strong3_bothfree,'black','o',color_map,norm,zorder=2)
plot_BvsDP_Kcolored(tmin,strong3_fixedSP,'black','v',color_map,norm,zorder=3)
plot_BvsDP_Kcolored(tmin,strong3_fixedOP,'black','^',color_map,norm,zorder=2)

# axis stuff
plt.xlim(-10,  37.5); plt.ylim(-10,  37.5)
plt.plot([-10, 37.5], [-10, 37.5], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30] )
fixed_aspect_ratio(1)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
# plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


