#!/bin/python
import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions_plotting import plot_BvsFullForce_wKthresh, plot_BvsDP_wKthresh, plot_BvsFullForce_wKthresh_overturned, plot_BvsDP_wKthresh_overturned
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
    plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.200-to-400kmB.no-OT.png'])
    plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.200-to-400kmB.no-OT.pdf'])
elif depths == "other":
    analysis_depth1 = 250e3
    analysis_depth2 = 300e3
    analysis_depth3 = 350e3
    plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.250-to-350kmB.no-OT.png'])
    plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.250-to-350kmB.no-OT.pdf'])
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
text1_ref_bothfree  = ''.join(['text_files/TESTB/',name_ref_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedSP   = ''.join(['text_files/TESTB/',name_ref_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedOP   = ''.join(['text_files/TESTB/',name_ref_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_bothfree  = ''.join(['text_files/TESTB/',name_ref2_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_fixedSP   = ''.join(['text_files/TESTB/',name_ref2_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_fixedOP   = ''.join(['text_files/TESTB/',name_ref2_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_bothfree = ''.join(['text_files/TESTB/',name_weak_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedSP      = ''.join(['text_files/TESTB/',name_weak_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedOP      = ''.join(['text_files/TESTB/',name_weak_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_bothfree= ''.join(['text_files/TESTB/',name_strong_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedSP = ''.join(['text_files/TESTB/',name_strong_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedOP    = ''.join(['text_files/TESTB/',name_strong_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_bothfree=''.join(['text_files/TESTB/',name_new_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_fixedSP= ''.join(['text_files/TESTB/',name_new_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_fixedOP= ''.join(['text_files/TESTB/',name_new_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# intermediate
text2_ref_bothfree  = ''.join(['text_files/TESTB/',name_ref_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedSP   = ''.join(['text_files/TESTB/',name_ref_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedOP   = ''.join(['text_files/TESTB/',name_ref_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_bothfree  = ''.join(['text_files/TESTB/',name_ref2_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_fixedSP   = ''.join(['text_files/TESTB/',name_ref2_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_fixedOP   = ''.join(['text_files/TESTB/',name_ref2_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_bothfree = ''.join(['text_files/TESTB/',name_weak_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedSP      = ''.join(['text_files/TESTB/',name_weak_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedOP      = ''.join(['text_files/TESTB/',name_weak_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_bothfree= ''.join(['text_files/TESTB/',name_strong_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedSP = ''.join(['text_files/TESTB/',name_strong_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedOP    = ''.join(['text_files/TESTB/',name_strong_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_bothfree=''.join(['text_files/TESTB/',name_new_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_fixedOP= ''.join(['text_files/TESTB/',name_new_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_fixedSP= ''.join(['text_files/TESTB/',name_new_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# deep
text3_ref_bothfree  = ''.join(['text_files/TESTB/',name_ref_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedSP   = ''.join(['text_files/TESTB/',name_ref_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedOP   = ''.join(['text_files/TESTB/',name_ref_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_bothfree  = ''.join(['text_files/TESTB/',name_ref2_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_fixedSP   = ''.join(['text_files/TESTB/',name_ref2_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_fixedOP   = ''.join(['text_files/TESTB/',name_ref2_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_bothfree = ''.join(['text_files/TESTB/',name_weak_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedSP      = ''.join(['text_files/TESTB/',name_weak_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedOP      = ''.join(['text_files/TESTB/',name_weak_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_bothfree= ''.join(['text_files/TESTB/',name_strong_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedSP = ''.join(['text_files/TESTB/',name_strong_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedOP    = ''.join(['text_files/TESTB/',name_strong_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_bothfree=''.join(['text_files/TESTB/',name_new_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_fixedSP= ''.join(['text_files/TESTB/',name_new_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_fixedOP= ''.join(['text_files/TESTB/',name_new_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


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


# DEPTH 1
ax=fig.add_subplot(gs[0,0])
# 50
plot_BvsFullForce_wKthresh(tmin,weak1_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak1_fixedSP,curvature_thresh, 'tan','black','v',zorder=4)
plot_BvsFullForce_wKthresh(tmin,weak1_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref1_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref1_fixedSP,curvature_thresh, 'peru','black','v',zorder=4)
plot_BvsFullForce_wKthresh(tmin,ref1_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new1_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new1_fixedSP,curvature_thresh, 'firebrick','black','v',zorder=4)
plot_BvsFullForce_wKthresh_overturned(tmin,new1_fixedOP,curvature_thresh, 'firebrick','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref1b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref1b_fixedSP,curvature_thresh, 'maroon','black','v',zorder=4)
plot_BvsFullForce_wKthresh_overturned(tmin,ref1b_fixedOP,curvature_thresh, 'maroon','^')
# 1000
plot_BvsFullForce_wKthresh_overturned(tmin,strong1_bothfree,curvature_thresh,'black','o')
plot_BvsFullForce_wKthresh(tmin,strong1_fixedSP,curvature_thresh, 'black','black','v',zorder=2)
plot_BvsFullForce_wKthresh_overturned(tmin,strong1_fixedOP,curvature_thresh, 'black','^')

# axis stuff
plt.xlim(-3,  30); plt.ylim(-3,  30)
plt.plot([-3, 30], [-3, 30], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([ 0, 10, 20, 30])
ax.set_xticklabels([ 0, 10, 20, 30])
ax.set_yticks([ 0, 10, 20, 30])
ax.set_yticklabels([ 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)


fixed_aspect_ratio(1)


# DEPTH 2
ax=fig.add_subplot(gs[0,1])
# 50
plot_BvsFullForce_wKthresh(tmin,weak2_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak2_fixedSP,curvature_thresh, 'tan','black','v',zorder=4)
plot_BvsFullForce_wKthresh(tmin,weak2_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref2_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref2_fixedSP,curvature_thresh, 'peru','black','v',zorder=4)
plot_BvsFullForce_wKthresh(tmin,ref2_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new2_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new2_fixedSP,curvature_thresh, 'firebrick','black','v',zorder=4)
plot_BvsFullForce_wKthresh_overturned(tmin,new2_fixedOP,curvature_thresh, 'firebrick','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref2b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref2b_fixedSP,curvature_thresh, 'maroon','black','v',zorder=4)
plot_BvsFullForce_wKthresh_overturned(tmin,ref2b_fixedOP,curvature_thresh, 'maroon','^')
# 1000
plot_BvsFullForce_wKthresh_overturned(tmin,strong2_bothfree,curvature_thresh,'black','o')
plot_BvsFullForce_wKthresh(tmin,strong2_fixedSP,curvature_thresh, 'black','black','v',zorder=2)
plot_BvsFullForce_wKthresh_overturned(tmin,strong2_fixedOP,curvature_thresh, 'black','^')

# axis stuff
plt.xlim(-3,  30); plt.ylim(-3,  30)
plt.plot([-3, 30], [-3, 30], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([ 0, 10, 20, 30])
ax.set_xticklabels([ 0, 10, 20, 30])
ax.set_yticks([ 0, 10, 20, 30])
ax.set_yticklabels([ 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)

fixed_aspect_ratio(1)



# DEPTH 3
ax=fig.add_subplot(gs[0,2])
# 50
plot_BvsFullForce_wKthresh(tmin,weak3_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak3_fixedSP,curvature_thresh, 'tan','black','v',zorder=4)
plot_BvsFullForce_wKthresh(tmin,weak3_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref3_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref3_fixedSP,curvature_thresh, 'peru','black','v',zorder=4)
plot_BvsFullForce_wKthresh(tmin,ref3_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new3_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new3_fixedSP,curvature_thresh, 'firebrick','black','v',zorder=4)
plot_BvsFullForce_wKthresh_overturned(tmin,new3_fixedOP,curvature_thresh, 'firebrick','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref3b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref3b_fixedSP,curvature_thresh, 'maroon','black','v',zorder=4)
plot_BvsFullForce_wKthresh_overturned(tmin,ref3b_fixedOP,curvature_thresh, 'maroon','^')
# 1000
plot_BvsFullForce_wKthresh_overturned(tmin,strong3_bothfree,curvature_thresh,'black','o')
plot_BvsFullForce_wKthresh(tmin,strong3_fixedSP,curvature_thresh, 'black','black','v',zorder=2)
plot_BvsFullForce_wKthresh_overturned(tmin,strong3_fixedOP,curvature_thresh, 'black','^')

# axis stuff
# axis stuff
plt.xlim(-3,  30); plt.ylim(-3,  30)
plt.plot([-3, 30], [-3, 30], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([ 0, 10, 20, 30])
ax.set_xticklabels([ 0, 10, 20, 30])
ax.set_yticks([ 0, 10, 20, 30])
ax.set_yticklabels([ 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)

fixed_aspect_ratio(1)


#########################################################################################
################################ Now just DP vs. DP #####################################
#########################################################################################

# DEPTH 1
ax=fig.add_subplot(gs[1,0])
# 50
plot_BvsDP_wKthresh(tmin,weak1_bothfree,curvature_thresh,'tan','black','o')
plot_BvsDP_wKthresh(tmin,weak1_fixedSP,curvature_thresh, 'tan','black','v',zorder=4)
plot_BvsDP_wKthresh(tmin,weak1_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsDP_wKthresh(tmin,ref1_bothfree,curvature_thresh,'peru','black','o')
plot_BvsDP_wKthresh(tmin,ref1_fixedSP,curvature_thresh, 'peru','black','v',zorder=4)
plot_BvsDP_wKthresh(tmin,ref1_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsDP_wKthresh(tmin,new1_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsDP_wKthresh(tmin,new1_fixedSP,curvature_thresh, 'firebrick','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,new1_fixedOP,curvature_thresh, 'firebrick','^')
# 500
plot_BvsDP_wKthresh(tmin,ref1b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsDP_wKthresh(tmin,ref1b_fixedSP,curvature_thresh, 'maroon','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,ref1b_fixedOP,curvature_thresh, 'maroon','^')
# 1000
plot_BvsDP_wKthresh_overturned(tmin,strong1_bothfree,curvature_thresh,'black','o')
plot_BvsDP_wKthresh(tmin,strong1_fixedSP,curvature_thresh, 'black','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,strong1_fixedOP,curvature_thresh, 'black','^')

# axis stuff
plt.xlim(-3,  30); plt.ylim(-3,  30)
plt.plot([-3, 30], [-3, 30], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([ 0, 10, 20, 30])
ax.set_xticklabels([ 0, 10, 20, 30])
ax.set_yticks([ 0, 10, 20, 30])
ax.set_yticklabels([ 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)

fixed_aspect_ratio(1)


# DEPTH 2
ax=fig.add_subplot(gs[1,1])
# 50
plot_BvsDP_wKthresh(tmin,weak2_bothfree,curvature_thresh,'tan','black','o')
plot_BvsDP_wKthresh(tmin,weak2_fixedSP,curvature_thresh, 'tan','black','v',zorder=4)
plot_BvsDP_wKthresh(tmin,weak2_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsDP_wKthresh(tmin,ref2_bothfree,curvature_thresh,'peru','black','o')
plot_BvsDP_wKthresh(tmin,ref2_fixedSP,curvature_thresh, 'peru','black','v',zorder=4)
plot_BvsDP_wKthresh(tmin,ref2_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsDP_wKthresh(tmin,new2_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsDP_wKthresh(tmin,new2_fixedSP,curvature_thresh, 'firebrick','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,new2_fixedOP,curvature_thresh, 'firebrick','^')
# 500
plot_BvsDP_wKthresh(tmin,ref2b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsDP_wKthresh(tmin,ref2b_fixedSP,curvature_thresh, 'maroon','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,ref2b_fixedOP,curvature_thresh, 'maroon','^')
# 1000
plot_BvsDP_wKthresh_overturned(tmin,strong2_bothfree,curvature_thresh,'black','o')
plot_BvsDP_wKthresh(tmin,strong2_fixedSP,curvature_thresh, 'black','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,strong2_fixedOP,curvature_thresh, 'black','^')

# axis stuff
plt.xlim(-3,  30); plt.ylim(-3,  30)
plt.plot([-3, 30], [-3, 30], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([ 0, 10, 20, 30])
ax.set_xticklabels([ 0, 10, 20, 30])
ax.set_yticks([ 0, 10, 20, 30])
ax.set_yticklabels([ 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)

fixed_aspect_ratio(1)

# DEPTH 3
ax=fig.add_subplot(gs[1,2])
# 50
plot_BvsDP_wKthresh(tmin,weak3_bothfree,curvature_thresh,'tan','black','o')
plot_BvsDP_wKthresh(tmin,weak3_fixedSP,curvature_thresh, 'tan','black','v',zorder=4)
plot_BvsDP_wKthresh(tmin,weak3_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsDP_wKthresh(tmin,ref3_bothfree,curvature_thresh,'peru','black','o')
plot_BvsDP_wKthresh(tmin,ref3_fixedSP,curvature_thresh, 'peru','black','v',zorder=4)
plot_BvsDP_wKthresh(tmin,ref3_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsDP_wKthresh(tmin,new3_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsDP_wKthresh(tmin,new3_fixedSP,curvature_thresh, 'firebrick','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,new3_fixedOP,curvature_thresh, 'firebrick','^')
# 500
plot_BvsDP_wKthresh(tmin,ref3b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsDP_wKthresh(tmin,ref3b_fixedSP,curvature_thresh, 'maroon','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,ref3b_fixedOP,curvature_thresh, 'maroon','^')
# 1000
plot_BvsDP_wKthresh_overturned(tmin,strong3_bothfree,curvature_thresh,'black','o')
plot_BvsDP_wKthresh(tmin,strong3_fixedSP,curvature_thresh, 'black','black','v',zorder=4)
plot_BvsDP_wKthresh_overturned(tmin,strong3_fixedOP,curvature_thresh, 'black','^')


# axis stuff
plt.xlim(-3,  30); plt.ylim(-3,  30)
plt.plot([-3, 30], [-3, 30], color='black', linewidth=1, zorder=2)
plt.ylabel(r'$\Delta P$   [MPa]')
plt.xlabel(r'$B_{slab}$   [MPa]')
ax.set_xticks([ 0, 10, 20, 30])
ax.set_xticklabels([ 0, 10, 20, 30])
ax.set_yticks([ 0, 10, 20, 30])
ax.set_yticklabels([ 0, 10, 20, 30])
plt.minorticks_on()
plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=1)

fixed_aspect_ratio(1)

plt.subplots_adjust(wspace=0.3, hspace=0.3)
plt.savefig(plot_name_png, format='png', dpi=600, bbox_inches='tight')
plt.savefig(plot_name_pdf, format='pdf', bbox_inches='tight')


