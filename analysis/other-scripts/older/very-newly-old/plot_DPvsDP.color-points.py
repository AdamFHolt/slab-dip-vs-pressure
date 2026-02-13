#!/bin/python
import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions_plotting import plot_BvsFullForce_wKthresh, plot_BvsDP_wKthresh
from functions_plotting import plot_BvsFullForce_Kcolored, plot_BvsDP_scalingcolored
import matplotlib.font_manager as fm
font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
myriad_pro = fm.FontProperties(fname=font_path)

mpl.rcParams['font.family'] = 'Myriad Pro'  # Now it should work if properly installed!
mpl.rcParams['font.size'] = 7.5
mpl.rcParams['axes.labelsize'] = 7.5
mpl.rcParams['axes.labelpad'] = 1.25
mpl.rcParams['xtick.labelsize'] = 6.5
mpl.rcParams['ytick.labelsize'] = 6.5
mpl.rcParams['xtick.major.pad'] = 2
mpl.rcParams['ytick.major.pad'] = 2
mpl.rcParams['xtick.major.size'] = 2.5
mpl.rcParams['ytick.major.size'] = 2.5
mpl.rcParams['xtick.minor.size'] = 1.25
mpl.rcParams['ytick.minor.size'] = 1.25

# print("Matplotlib is using:", plt.rcParams['font.family'])

analysis_depth_dz = float(sys.argv[1])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[2])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[3])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)
depths = str(sys.argv[4])               # "normal" (230e3, 330e3, 430e3) or "other" (250e3, 300e3, 350e3)
curvature_thresh = 1e9

tactual_min = 11 # first time step to use
tmin = tactual_min - 8

if depths == "normal":
    analysis_depth1 = 200e3
    analysis_depth2 = 300e3
    analysis_depth3 = 400e3
    plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.200-to-400km.color-points.png'])
    plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.200-to-400km.color-points.pdf'])
elif depths == "other":
    analysis_depth1 = 250e3
    analysis_depth2 = 300e3
    analysis_depth3 = 350e3
    plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.250-to-350km.color-points.png'])
    plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.250-to-350km.color-points.pdf'])
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
new3_bothfree = np.loadtxt((text3_new_bothfree)) 
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
	im = plt.scatter(weak1_bothfree[i,4]/1.e6,(weak1_bothfree[i,3]-weak1_bothfree[i,6]+weak1_bothfree[i,17])/1.e6,s=10,c=np.abs(weak1_bothfree[i,11])*1e3,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.1,zorder=3,marker='o')
plot_BvsFullForce_Kcolored(tmin,weak1_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,weak1_fixedOP,'black','^',color_map,norm,zorder=3)
# 250
plot_BvsFullForce_Kcolored(tmin,ref1_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref1_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,ref1_fixedOP,'black','^',color_map,norm,zorder=3)
# 375
plot_BvsFullForce_Kcolored(tmin,new1_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,new1_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,new1_fixedOP,'black','^',color_map,norm,zorder=3)
# 500
plot_BvsFullForce_Kcolored(tmin,ref1b_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref1b_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,ref1b_fixedOP,'black','^',color_map,norm,zorder=3)
# 1000
plot_BvsFullForce_Kcolored(tmin,strong1_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,strong1_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,strong1_fixedOP,'black','^',color_map,norm,zorder=3)


# axis stuff
plt.xlim(-5,  34); plt.ylim(-5,  34)
plt.plot([-5, 34], [-5, 34], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([-5, 0, 10, 20, 30])
ax.set_xticklabels([-5, 0, 10, 20, 30])
ax.set_yticks([-5, 0, 10, 20, 30])
ax.set_yticklabels([-5, 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)

# c-bar
cbar = plt.colorbar(im, cax = fig.add_axes([0.31, 0.585, 0.008, 0.08]), orientation='vertical', ticks=[0,0.001,0.002], extend='max') # left, base, length left-right, length top-bott
cbar.ax.tick_params(axis='y',labelsize=4.5,pad=1)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.yaxis.set_label_position('left')
cbar.set_label(r'$K$  [1/km]', rotation=0, labelpad=-10, y=1.1)
cbar.ax.tick_params(axis='y', labelsize=5.5, pad=1)

# DEPTH 2
ax=fig.add_subplot(gs[0,1])
# 50
plot_BvsFullForce_Kcolored(tmin,weak2_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,weak2_fixedSP,'black','v',color_map,norm, zorder=4)
plot_BvsFullForce_Kcolored(tmin,weak2_fixedOP,'black','^',color_map,norm, zorder=3)
# 250
plot_BvsFullForce_Kcolored(tmin,ref2_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref2_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,ref2_fixedOP,'black','^',color_map,norm,zorder=3)
# 375
plot_BvsFullForce_Kcolored(tmin,new2_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,new2_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,new2_fixedOP,'black','^',color_map,norm,zorder=3)
# 500
plot_BvsFullForce_Kcolored(tmin,ref2b_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref2b_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,ref2b_fixedOP,'black','^',color_map,norm,zorder=3)
# 1000
plot_BvsFullForce_Kcolored(tmin,strong2_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,strong2_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,strong2_fixedOP,'black','^',color_map,norm,zorder=3)

# axis stuff
plt.xlim(-5,  34); plt.ylim(-5,  34)
plt.plot([-5, 34], [-5, 34], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([-5, 0, 10, 20, 30])
ax.set_xticklabels([-5, 0, 10, 20, 30])
ax.set_yticks([-5, 0, 10, 20, 30])
ax.set_yticklabels([-5, 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)

# DEPTH 3
ax=fig.add_subplot(gs[0,2])
# 50
plot_BvsFullForce_Kcolored(tmin,weak3_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,weak3_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,weak3_fixedOP,'black','^',color_map,norm,zorder=3)
# 250
plot_BvsFullForce_Kcolored(tmin,ref3_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref3_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,ref3_fixedOP,'black','^',color_map,norm,zorder=3)
# 375
plot_BvsFullForce_Kcolored(tmin,new3_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,new3_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,new3_fixedOP,'black','^',color_map,norm,zorder=3)
# 500
plot_BvsFullForce_Kcolored(tmin,ref3b_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,ref3b_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,ref3b_fixedOP,'black','^',color_map,norm,zorder=3)
# 1000
plot_BvsFullForce_Kcolored(tmin,strong3_bothfree,'black','o',color_map,norm,zorder=3)
plot_BvsFullForce_Kcolored(tmin,strong3_fixedSP,'black','v',color_map,norm,zorder=4)
plot_BvsFullForce_Kcolored(tmin,strong3_fixedOP,'black','^',color_map,norm,zorder=3)

# axis stuff
# axis stuff
plt.xlim(-5,  34); plt.ylim(-5,  34)
plt.plot([-5, 34], [-5, 34], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([-5, 0, 10, 20, 30])
ax.set_xticklabels([-5, 0, 10, 20, 30])
ax.set_yticks([-5, 0, 10, 20, 30])
ax.set_yticklabels([-5, 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)


#########################################################################################
################################ Now just DP vs. DP #####################################
#########################################################################################



# set up color map
color_map2 = cm.get_cmap('viridis')
color_min2 = 0
color_max2 = 200
interval2  = color_max2/14.
norm2 = matplotlib.colors.BoundaryNorm(np.arange(color_min2,color_max2+interval2,interval2), color_map.N)
mant_visc = 2.5e20

# DEPTH 1
ax=fig.add_subplot(gs[1,0])

# 50
cmyr_to_ms = 0.01/(365.25*24*3600)
for i in range(tmin,len(weak1_bothfree)):
	im2 = plt.scatter(weak1_bothfree[i,4]/1.e6,weak1_bothfree[i,3]/1.e6,s=10,c=weak1_bothfree[i,11]*50*mant_visc*weak1_bothfree[i,19]*cmyr_to_ms*1e-6,cmap=color_map2,norm=norm2,edgecolor='black',linewidth=0.1,zorder=3,marker='o')
plot_BvsDP_scalingcolored(tmin,weak1_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=50)
plot_BvsDP_scalingcolored(tmin,weak1_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=50)
# 250
plot_BvsDP_scalingcolored(tmin,ref1_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
plot_BvsDP_scalingcolored(tmin,ref1_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=250)
plot_BvsDP_scalingcolored(tmin,ref1_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
# 375
plot_BvsDP_scalingcolored(tmin,new1_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=375)
plot_BvsDP_scalingcolored(tmin,new1_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=375)
plot_BvsDP_scalingcolored(tmin,new1_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=375)
# 500
plot_BvsDP_scalingcolored(tmin,ref1b_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=500)
plot_BvsDP_scalingcolored(tmin,ref1b_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=500)
plot_BvsDP_scalingcolored(tmin,ref1b_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=500)
# 1000
plot_BvsDP_scalingcolored(tmin,strong1_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=1000)
plot_BvsDP_scalingcolored(tmin,strong1_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=1000)
plot_BvsDP_scalingcolored(tmin,strong1_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=1000)



# axis stuff
plt.xlim(-5,  34); plt.ylim(-5,  34)
plt.plot([-5, 34], [-5, 34], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([-5, 0, 10, 20, 30])
ax.set_xticklabels([-5, 0, 10, 20, 30])
ax.set_yticks([-5, 0, 10, 20, 30])
ax.set_yticklabels([-5, 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)


fixed_aspect_ratio(1)

# c-bar
cbar = plt.colorbar(im2, cax = fig.add_axes([0.18, 0.272, 0.008, 0.08]), orientation='vertical', ticks=[0,color_max2/2,color_max2], extend='max') # left, base, length left-right, length top-bott
cbar.ax.tick_params(axis='y',labelsize=4.5,pad=1)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.yaxis.set_label_position('left')
cbar.set_label(r'$\eta K V_{C}$  [MPa]', rotation=0, labelpad=-10, y=1.1)
cbar.ax.tick_params(axis='y', labelsize=5.5, pad=1)


# DEPTH 2
ax=fig.add_subplot(gs[1,1])
# 50
plot_BvsDP_scalingcolored(tmin,weak2_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
plot_BvsDP_scalingcolored(tmin,weak2_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=50)
plot_BvsDP_scalingcolored(tmin,weak2_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=50)
# 250
plot_BvsDP_scalingcolored(tmin,ref2_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
plot_BvsDP_scalingcolored(tmin,ref2_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=250)
plot_BvsDP_scalingcolored(tmin,ref2_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
# 375
plot_BvsDP_scalingcolored(tmin,new2_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=375)
plot_BvsDP_scalingcolored(tmin,new2_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=375)
plot_BvsDP_scalingcolored(tmin,new2_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=375)
# 500
plot_BvsDP_scalingcolored(tmin,ref2b_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=500)
plot_BvsDP_scalingcolored(tmin,ref2b_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=500)
plot_BvsDP_scalingcolored(tmin,ref2b_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=500)
# 1000
plot_BvsDP_scalingcolored(tmin,strong2_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=1000)
plot_BvsDP_scalingcolored(tmin,strong2_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=1000)
plot_BvsDP_scalingcolored(tmin,strong2_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=1000)


# axis stuff
plt.xlim(-5,  34); plt.ylim(-5,  34)
plt.plot([-5, 34], [-5, 34], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([-5, 0, 10, 20, 30])
ax.set_xticklabels([-5, 0, 10, 20, 30])
ax.set_yticks([-5, 0, 10, 20, 30])
ax.set_yticklabels([-5, 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)

fixed_aspect_ratio(1)

# DEPTH 3
ax=fig.add_subplot(gs[1,2])
# 50
plot_BvsDP_scalingcolored(tmin,weak3_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
plot_BvsDP_scalingcolored(tmin,weak3_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=50)
plot_BvsDP_scalingcolored(tmin,weak3_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=50)
# 250
plot_BvsDP_scalingcolored(tmin,ref3_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
plot_BvsDP_scalingcolored(tmin,ref3_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=250)
plot_BvsDP_scalingcolored(tmin,ref3_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=250)
# 375
plot_BvsDP_scalingcolored(tmin,new3_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=375)
plot_BvsDP_scalingcolored(tmin,new3_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=375)
plot_BvsDP_scalingcolored(tmin,new3_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=375)
# 500
plot_BvsDP_scalingcolored(tmin,ref3b_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=500)
plot_BvsDP_scalingcolored(tmin,ref3b_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=500)
plot_BvsDP_scalingcolored(tmin,ref3b_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=500)
# 1000
plot_BvsDP_scalingcolored(tmin,strong3_bothfree,'black','o',color_map2,norm2,mant_visc,zorder=3,viscosity=1000)
plot_BvsDP_scalingcolored(tmin,strong3_fixedSP,'black','v',color_map2,norm2,mant_visc,zorder=4,viscosity=1000)
plot_BvsDP_scalingcolored(tmin,strong3_fixedOP,'black','^',color_map2,norm2,mant_visc,zorder=3,viscosity=1000)


# axis stuff
plt.xlim(-5,  34); plt.ylim(-5,  34)
plt.plot([-5, 34], [-5, 34], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([-5, 0, 10, 20, 30])
ax.set_xticklabels([-5, 0, 10, 20, 30])
ax.set_yticks([-5, 0, 10, 20, 30])
ax.set_yticklabels([-5, 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)

fixed_aspect_ratio(1)

plt.subplots_adjust(wspace=0.3, hspace=0.3)
plt.savefig(plot_name_png, format='png', dpi=600, bbox_inches='tight')
plt.savefig(plot_name_pdf, format='pdf', bbox_inches='tight')


