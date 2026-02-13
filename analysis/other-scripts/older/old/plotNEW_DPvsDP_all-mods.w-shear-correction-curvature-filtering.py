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

analysis_depth_dz = float(sys.argv[1])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[2])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[3])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

analysis_depth1 = 230e3
analysis_depth2 = 330e3
analysis_depth3 = 430e3


name_weak_bothfree   = "2D_compositional_subd_lower-res_new_50plates"
name_weak_fixedSP    = "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name_weak_fixedOP    = "2D_compositional_subd_lower-res_new_FixedOP_50plates"
name_ref_bothfree    = "2D_compositional_subd_lower-res_new_250plates"
name_ref_fixedSP     = "2D_compositional_subd_lower-res_new_FixedSP_250plates"
name_ref_fixedOP     = "2D_compositional_subd_lower-res_new_FixedOP_250plates"
name_ref2_bothfree   = "2D_compositional_subd_lower-res_new"         # 500
name_ref2_fixedSP    = "2D_compositional_subd_FixedSP_lower-res_new" # 500 
name_ref2_fixedOP    = "2D_compositional_subd_FixedOP_lower-res_new" # 500
name_strong_bothfree = "2D_compositional_subd_lower-res_new_1000plates"
name_strong_fixedSP  = "2D_compositional_subd_lower-res_new_FixedSP_1000plates"
name_strong_fixedOP  = "2D_compositional_subd_lower-res_new_FixedOP_1000plates"
name_new_bothfree    = "2D_compositional_subd_lower-res_new_375plates"
name_new_fixedSP     = "2D_compositional_subd_lower-res_new_FixedSP_375plates"
name_new_fixedOP     = "2D_compositional_subd_lower-res_new_FixedOP_375plates"

# shall depths
text1_ref_bothfree  = ''.join(['text_files/new/',name_ref_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedSP   = ''.join(['text_files/new/',name_ref_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedOP   = ''.join(['text_files/new/',name_ref_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_bothfree  = ''.join(['text_files/new/',name_ref2_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_fixedSP   = ''.join(['text_files/new/',name_ref2_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref2_fixedOP   = ''.join(['text_files/new/',name_ref2_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_bothfree = ''.join(['text_files/new/',name_weak_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedSP      = ''.join(['text_files/new/',name_weak_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedOP      = ''.join(['text_files/new/',name_weak_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_bothfree= ''.join(['text_files/new/',name_strong_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedSP = ''.join(['text_files/new/',name_strong_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedOP    = ''.join(['text_files/new/',name_strong_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_bothfree=''.join(['text_files/new/',name_new_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_fixedSP= ''.join(['text_files/new/',name_new_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_new_fixedOP= ''.join(['text_files/new/',name_new_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# intermediate 
text2_ref_bothfree  = ''.join(['text_files/new/',name_ref_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedSP   = ''.join(['text_files/new/',name_ref_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedOP   = ''.join(['text_files/new/',name_ref_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_bothfree  = ''.join(['text_files/new/',name_ref2_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_fixedSP   = ''.join(['text_files/new/',name_ref2_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref2_fixedOP   = ''.join(['text_files/new/',name_ref2_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_bothfree = ''.join(['text_files/new/',name_weak_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedSP      = ''.join(['text_files/new/',name_weak_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedOP      = ''.join(['text_files/new/',name_weak_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_bothfree= ''.join(['text_files/new/',name_strong_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedSP = ''.join(['text_files/new/',name_strong_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedOP    = ''.join(['text_files/new/',name_strong_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_bothfree=''.join(['text_files/new/',name_new_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_fixedSP= ''.join(['text_files/new/',name_new_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_new_fixedOP= ''.join(['text_files/new/',name_new_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# deep
text3_ref_bothfree  = ''.join(['text_files/new/',name_ref_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedSP   = ''.join(['text_files/new/',name_ref_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedOP   = ''.join(['text_files/new/',name_ref_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_bothfree  = ''.join(['text_files/new/',name_ref2_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_fixedSP   = ''.join(['text_files/new/',name_ref2_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref2_fixedOP   = ''.join(['text_files/new/',name_ref2_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_bothfree = ''.join(['text_files/new/',name_weak_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedSP      = ''.join(['text_files/new/',name_weak_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedOP      = ''.join(['text_files/new/',name_weak_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_bothfree= ''.join(['text_files/new/',name_strong_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedSP = ''.join(['text_files/new/',name_strong_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedOP    = ''.join(['text_files/new/',name_strong_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_bothfree=''.join(['text_files/new/',name_new_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_fixedSP= ''.join(['text_files/new/',name_new_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_new_fixedOP= ''.join(['text_files/new/',name_new_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


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

plot_name_png = ''.join(['plots/DP-comparisons/compilations/NEW/DP-vs-DP_all-mods.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-shear-correction-curvature-filtering.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/NEW/DP-vs-DP_all-mods.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-shear-correction-curvature-filtering.pdf'])

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


ax=fig.add_subplot(gs[0,0])
curvature_thresh = 0.00075 # 1/km

# 50
for i in range(3,len(weak1_bothfree)):
    if np.abs(weak1_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak1_bothfree[i,4]/1.e6,(weak1_bothfree[i,3]-weak1_bothfree[i,6])/1.e6,color='tan',s=10,edgecolor='black',lw=0.1,zorder=3) 
    else:
        plt.scatter(weak1_bothfree[i,4]/1.e6,(weak1_bothfree[i,3]-weak1_bothfree[i,6])/1.e6,color='tan',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(weak1_fixedSP)):
    if np.abs(weak1_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak1_fixedSP[i,4]/1.e6,(weak1_fixedSP[i,3]-weak1_fixedSP[i,6])/1.e6,color='tan',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(weak1_fixedSP[i,4]/1.e6,(weak1_fixedSP[i,3]-weak1_fixedSP[i,6])/1.e6,color='tan',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(weak1_fixedOP)):
    if np.abs(weak1_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak1_fixedOP[i,4]/1.e6,(weak1_fixedOP[i,3]-weak1_fixedOP[i,6])/1.e6,color='tan',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3) 
    else:
        plt.scatter(weak1_fixedOP[i,4]/1.e6,(weak1_fixedOP[i,3]-weak1_fixedOP[i,6])/1.e6,color='tan',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 250
for i in range(3,len(ref1_bothfree)):
    if np.abs(ref1_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref1_bothfree[i,4]/1.e6,(ref1_bothfree[i,3]-ref1_bothfree[i,6])/1.e6,color='peru',s=10,edgecolor='black',lw=0.1,zorder=3) 
    else:
        plt.scatter(ref1_bothfree[i,4]/1.e6,(ref1_bothfree[i,3]-ref1_bothfree[i,6])/1.e6,color='peru',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref1_fixedSP)):
    if np.abs(ref1_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref1_fixedSP[i,4]/1.e6,(ref1_fixedSP[i,3]-ref1_fixedSP[i,6])/1.e6,color='peru',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3) 
    else:
        plt.scatter(ref1_fixedSP[i,4]/1.e6,(ref1_fixedSP[i,3]-ref1_fixedSP[i,6])/1.e6,color='peru',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref1_fixedOP)):
    if np.abs(ref1_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref1_fixedOP[i,4]/1.e6,(ref1_fixedOP[i,3]-ref1_fixedOP[i,6])/1.e6,color='peru',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else:
        plt.scatter(ref1_fixedOP[i,4]/1.e6,(ref1_fixedOP[i,3]-ref1_fixedOP[i,6])/1.e6,color='peru',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 
# 375
for i in range(3,len(new1_bothfree)):
    if np.abs(new1_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(new1_bothfree[i,4]/1.e6,(new1_bothfree[i,3]-new1_bothfree[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',lw=0.1,zorder=3) 
    else:
        plt.scatter(new1_bothfree[i,4]/1.e6,(new1_bothfree[i,3]-new1_bothfree[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(new1_fixedSP)):
    if np.abs(new1_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(new1_fixedSP[i,4]/1.e6,(new1_fixedSP[i,3]-new1_fixedSP[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3) 
    else:
        plt.scatter(new1_fixedSP[i,4]/1.e6,(new1_fixedSP[i,3]-new1_fixedSP[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(new1_fixedOP)):
    if np.abs(new1_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(new1_fixedOP[i,4]/1.e6,(new1_fixedOP[i,3]-new1_fixedOP[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else:
        plt.scatter(new1_fixedOP[i,4]/1.e6,(new1_fixedOP[i,3]-new1_fixedOP[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 
# 500
for i in range(3,len(ref1b_bothfree)):
    if np.abs(ref1b_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref1b_bothfree[i,4]/1.e6,(ref1b_bothfree[i,3]-ref1b_bothfree[i,6])/1.e6,color='maroon',s=10,edgecolor='black',lw=0.1,zorder=3) 
    else:
        plt.scatter(ref1b_bothfree[i,4]/1.e6,(ref1b_bothfree[i,3]-ref1b_bothfree[i,6])/1.e6,color='maroon',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref1b_fixedSP)):
    if np.abs(ref1b_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref1b_fixedSP[i,4]/1.e6,(ref1b_fixedSP[i,3]-ref1b_fixedSP[i,6])/1.e6,color='maroon',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3) 
    else:
        plt.scatter(ref1b_fixedSP[i,4]/1.e6,(ref1b_fixedSP[i,3]-ref1b_fixedSP[i,6])/1.e6,color='maroon',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref1b_fixedOP)):
    if np.abs(ref1b_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref1b_fixedOP[i,4]/1.e6,(ref1b_fixedOP[i,3]-ref1b_fixedOP[i,6])/1.e6,color='maroon',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else:
        plt.scatter(ref1b_fixedOP[i,4]/1.e6,(ref1b_fixedOP[i,3]-ref1b_fixedOP[i,6])/1.e6,color='maroon',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 
# 1000
for i in range(3,len(strong1_bothfree)):
    if np.abs(strong1_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong1_bothfree[i,4]/1.e6,(strong1_bothfree[i,3]-strong1_bothfree[i,6])/1.e6,color='black',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(strong1_bothfree[i,4]/1.e6,(strong1_bothfree[i,3]-strong1_bothfree[i,6])/1.e6,color='black',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(strong1_fixedSP)):
    if np.abs(strong1_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong1_fixedSP[i,4]/1.e6,(strong1_fixedSP[i,3]-strong1_fixedSP[i,6])/1.e6,color='black',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(strong1_fixedSP[i,4]/1.e6,(strong1_fixedSP[i,3]-strong1_fixedSP[i,6])/1.e6,color='black',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(strong1_fixedOP)):
    if np.abs(strong1_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong1_fixedOP[i,4]/1.e6,(strong1_fixedOP[i,3]-strong1_fixedOP[i,6])/1.e6,color='black',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3) 
    else:
        plt.scatter(strong1_fixedOP[i,4]/1.e6,(strong1_fixedOP[i,3]-strong1_fixedOP[i,6])/1.e6,color='black',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# axis stuff
plt.xlim(-15,  45); plt.ylim(-15,  45)
plt.plot([-15, 45], [-15, 45], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel(r'$\mathregular{\Delta P +}$ slab stress  [MPa]',size=6.5)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30,40] )
# plt.annotate('230 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=6.5,color='k')
fixed_aspect_ratio(1)


# misfit vs. curvature
ax=fig.add_subplot(gs[0,1])

# 50
for i in range(3,len(weak2_bothfree)):
    if np.abs(weak2_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak2_bothfree[i,4]/1.e6,(weak2_bothfree[i,3]-weak2_bothfree[i,6])/1.e6,color='tan',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(weak2_bothfree[i,4]/1.e6,(weak2_bothfree[i,3]-weak2_bothfree[i,6])/1.e6,color='tan',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(weak2_fixedSP)):
    if np.abs(weak2_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak2_fixedSP[i,4]/1.e6,(weak2_fixedSP[i,3]-weak2_fixedSP[i,6])/1.e6,color='tan',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3) 
    else:
        plt.scatter(weak2_fixedSP[i,4]/1.e6,(weak2_fixedSP[i,3]-weak2_fixedSP[i,6])/1.e6,color='tan',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(weak2_fixedOP)):
    if np.abs(weak2_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak2_fixedOP[i,4]/1.e6,(weak2_fixedOP[i,3]-weak2_fixedOP[i,6])/1.e6,color='tan',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3) 
    else:
        plt.scatter(weak2_fixedOP[i,4]/1.e6,(weak2_fixedOP[i,3]-weak2_fixedOP[i,6])/1.e6,color='tan',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 250
for i in range(3,len(ref2_bothfree)):
    if np.abs(ref2_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref2_bothfree[i,4]/1.e6,(ref2_bothfree[i,3]-ref2_bothfree[i,6])/1.e6,color='peru',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(ref2_bothfree[i,4]/1.e6,(ref2_bothfree[i,3]-ref2_bothfree[i,6])/1.e6,color='peru',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref2_fixedSP)):
    if np.abs(ref2_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref2_fixedSP[i,4]/1.e6,(ref2_fixedSP[i,3]-ref2_fixedSP[i,6])/1.e6,color='peru',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(ref2_fixedSP[i,4]/1.e6,(ref2_fixedSP[i,3]-ref2_fixedSP[i,6])/1.e6,color='peru',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref2_fixedOP)):
    if np.abs(ref2_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref2_fixedOP[i,4]/1.e6,(ref2_fixedOP[i,3]-ref2_fixedOP[i,6])/1.e6,color='peru',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3) 
    else:
        plt.scatter(ref2_fixedOP[i,4]/1.e6,(ref2_fixedOP[i,3]-ref2_fixedOP[i,6])/1.e6,color='peru',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 375
for i in range(3,len(new2_bothfree)):
    if np.abs(new2_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(new2_bothfree[i,4]/1.e6,(new2_bothfree[i,3]-new2_bothfree[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(new2_bothfree[i,4]/1.e6,(new2_bothfree[i,3]-new2_bothfree[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(new2_fixedSP)):
    if np.abs(new2_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(new2_fixedSP[i,4]/1.e6,(new2_fixedSP[i,3]-new2_fixedSP[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(new2_fixedSP[i,4]/1.e6,(new2_fixedSP[i,3]-new2_fixedSP[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(new2_fixedOP)):
    if np.abs(new2_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(new2_fixedOP[i,4]/1.e6,(new2_fixedOP[i,3]-new2_fixedOP[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3) 
    else:
        plt.scatter(new2_fixedOP[i,4]/1.e6,(new2_fixedOP[i,3]-new2_fixedOP[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 500
for i in range(3,len(ref2b_bothfree)):
    if np.abs(ref2b_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref2b_bothfree[i,4]/1.e6,(ref2b_bothfree[i,3]-ref2b_bothfree[i,6])/1.e6,color='maroon',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(ref2b_bothfree[i,4]/1.e6,(ref2b_bothfree[i,3]-ref2b_bothfree[i,6])/1.e6,color='maroon',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref2b_fixedSP)):
    if np.abs(ref2b_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref2b_fixedSP[i,4]/1.e6,(ref2b_fixedSP[i,3]-ref2b_fixedSP[i,6])/1.e6,color='maroon',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(ref2b_fixedSP[i,4]/1.e6,(ref2b_fixedSP[i,3]-ref2b_fixedSP[i,6])/1.e6,color='maroon',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref2b_fixedOP)):
    if np.abs(ref2b_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref2b_fixedOP[i,4]/1.e6,(ref2b_fixedOP[i,3]-ref2b_fixedOP[i,6])/1.e6,color='maroon',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3) 
    else:
        plt.scatter(ref2b_fixedOP[i,4]/1.e6,(ref2b_fixedOP[i,3]-ref2b_fixedOP[i,6])/1.e6,color='maroon',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 1000
for i in range(3,len(strong2_bothfree)):
    if np.abs(strong2_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong2_bothfree[i,4]/1.e6,(strong2_bothfree[i,3]-strong2_bothfree[i,6])/1.e6,color='black',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(strong2_bothfree[i,4]/1.e6,(strong2_bothfree[i,3]-strong2_bothfree[i,6])/1.e6,color='black',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(strong2_fixedSP)):
    if np.abs(strong2_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong2_fixedSP[i,4]/1.e6,(strong2_fixedSP[i,3]-strong2_fixedSP[i,6])/1.e6,color='black',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else: 
        plt.scatter(strong2_fixedSP[i,4]/1.e6,(strong2_fixedSP[i,3]-strong2_fixedSP[i,6])/1.e6,color='black',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(strong2_fixedOP)):
    if np.abs(strong2_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong2_fixedOP[i,4]/1.e6,(strong2_fixedOP[i,3]-strong2_fixedOP[i,6])/1.e6,color='black',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else: 
        plt.scatter(strong2_fixedOP[i,4]/1.e6,(strong2_fixedOP[i,3]-strong2_fixedOP[i,6])/1.e6,color='black',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# axis stuff
plt.xlim(-15,  45); plt.ylim(-15,  45)
plt.plot([-15, 45], [-15, 45], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30,40] )
# plt.annotate('330 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=6.5,color='k')
fixed_aspect_ratio(1)

# misfit vs. curvature
ax=fig.add_subplot(gs[0,2])

# 50
for i in range(3,len(weak3_bothfree)):
    if np.abs(weak3_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak3_bothfree[i,4]/1.e6,(weak3_bothfree[i,3]-weak3_bothfree[i,6])/1.e6,color='tan',s=10,edgecolor='black',lw=0.1,zorder=3)
    else: 
        plt.scatter(weak3_bothfree[i,4]/1.e6,(weak3_bothfree[i,3]-weak3_bothfree[i,6])/1.e6,color='tan',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(weak3_fixedSP)):
    if np.abs(weak3_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak3_fixedSP[i,4]/1.e6,(weak3_fixedSP[i,3]-weak3_fixedSP[i,6])/1.e6,color='tan',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else: 
        plt.scatter(weak3_fixedSP[i,4]/1.e6,(weak3_fixedSP[i,3]-weak3_fixedSP[i,6])/1.e6,color='tan',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(weak3_fixedOP)):
    if np.abs(weak3_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(weak3_fixedOP[i,4]/1.e6,(weak3_fixedOP[i,3]-weak3_fixedOP[i,6])/1.e6,color='tan',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else:
        plt.scatter(weak3_fixedOP[i,4]/1.e6,(weak3_fixedOP[i,3]-weak3_fixedOP[i,6])/1.e6,color='tan',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 250
for i in range(3,len(ref3_bothfree)):
    if np.abs(ref3_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref3_bothfree[i,4]/1.e6,(ref3_bothfree[i,3]-ref3_bothfree[i,6])/1.e6,color='peru',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(ref3_bothfree[i,4]/1.e6,(ref3_bothfree[i,3]-ref3_bothfree[i,6])/1.e6,color='peru',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref3_fixedSP)):
    if np.abs(ref3_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref3_fixedSP[i,4]/1.e6,(ref3_fixedSP[i,3]-ref3_fixedSP[i,6])/1.e6,color='peru',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(ref3_fixedSP[i,4]/1.e6,(ref3_fixedSP[i,3]-ref3_fixedSP[i,6])/1.e6,color='peru',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref3_fixedOP)):
    if np.abs(ref3_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref3_fixedOP[i,4]/1.e6,(ref3_fixedOP[i,3]-ref3_fixedOP[i,6])/1.e6,color='peru',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else: 
        plt.scatter(ref3_fixedOP[i,4]/1.e6,(ref3_fixedOP[i,3]-ref3_fixedOP[i,6])/1.e6,color='peru',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 375
for i in range(3,len(new3_bothfree)):
    if np.abs(new3_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(new3_bothfree[i,4]/1.e6,(new3_bothfree[i,3]-new3_bothfree[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(new3_bothfree[i,4]/1.e6,(new3_bothfree[i,3]-new3_bothfree[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(new3_fixedSP)):
    if np.abs(new3_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(new3_fixedSP[i,4]/1.e6,(new3_fixedSP[i,3]-new3_fixedSP[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(new3_fixedSP[i,4]/1.e6,(new3_fixedSP[i,3]-new3_fixedSP[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(new3_fixedOP)):
    if np.abs(new3_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(new3_fixedOP[i,4]/1.e6,(new3_fixedOP[i,3]-new3_fixedOP[i,6])/1.e6,color='firebrick',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else: 
        plt.scatter(new3_fixedOP[i,4]/1.e6,(new3_fixedOP[i,3]-new3_fixedOP[i,6])/1.e6,color='firebrick',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 500
for i in range(3,len(ref3b_bothfree)):
    if np.abs(ref3b_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref3b_bothfree[i,4]/1.e6,(ref3b_bothfree[i,3]-ref3b_bothfree[i,6])/1.e6,color='maroon',s=10,edgecolor='black',lw=0.1,zorder=3)
    else:
        plt.scatter(ref3b_bothfree[i,4]/1.e6,(ref3b_bothfree[i,3]-ref3b_bothfree[i,6])/1.e6,color='maroon',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref3b_fixedSP)):
    if np.abs(ref3b_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref3b_fixedSP[i,4]/1.e6,(ref3b_fixedSP[i,3]-ref3b_fixedSP[i,6])/1.e6,color='maroon',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else:
        plt.scatter(ref3b_fixedSP[i,4]/1.e6,(ref3b_fixedSP[i,3]-ref3b_fixedSP[i,6])/1.e6,color='maroon',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(ref3b_fixedOP)):
    if np.abs(ref3b_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(ref3b_fixedOP[i,4]/1.e6,(ref3b_fixedOP[i,3]-ref3b_fixedOP[i,6])/1.e6,color='maroon',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else: 
        plt.scatter(ref3b_fixedOP[i,4]/1.e6,(ref3b_fixedOP[i,3]-ref3b_fixedOP[i,6])/1.e6,color='maroon',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# 1000
for i in range(3,len(strong3_bothfree)):
    if np.abs(strong3_bothfree[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong3_bothfree[i,4]/1.e6,(strong3_bothfree[i,3]-strong3_bothfree[i,6])/1.e6,color='black',s=10,edgecolor='black',lw=0.1,zorder=3)
    else: 
        plt.scatter(strong3_bothfree[i,4]/1.e6,(strong3_bothfree[i,3]-strong3_bothfree[i,6])/1.e6,color='black',s=5,edgecolor='none',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(strong3_fixedSP)):
    if np.abs(strong3_fixedSP[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong3_fixedSP[i,4]/1.e6,(strong3_fixedSP[i,3]-strong3_fixedSP[i,6])/1.e6,color='black',s=10,edgecolor='black',marker='v',lw=0.1,zorder=3)
    else: 
        plt.scatter(strong3_fixedSP[i,4]/1.e6,(strong3_fixedSP[i,3]-strong3_fixedSP[i,6])/1.e6,color='black',s=5,edgecolor='none',marker='v',lw=0.1,zorder=2,alpha=0.2) 
for i in range(3,len(strong3_fixedOP)):
    if np.abs(strong3_fixedOP[i,11])*1000 < curvature_thresh: 
        plt.scatter(strong3_fixedOP[i,4]/1.e6,(strong3_fixedOP[i,3]-strong3_fixedOP[i,6])/1.e6,color='black',s=10,edgecolor='black',marker='^',lw=0.1,zorder=3)
    else: 
        plt.scatter(strong3_fixedOP[i,4]/1.e6,(strong3_fixedOP[i,3]-strong3_fixedOP[i,6])/1.e6,color='black',s=5,edgecolor='none',marker='^',lw=0.1,zorder=2,alpha=0.2) 

# axis stuff
plt.xlim(-15,  45); plt.ylim(-15,  45)
plt.plot([-15, 45], [-15, 45], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [-10,0,10,20,30,40] )
# plt.annotate('430 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=6.5,color='k')
fixed_aspect_ratio(1)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


