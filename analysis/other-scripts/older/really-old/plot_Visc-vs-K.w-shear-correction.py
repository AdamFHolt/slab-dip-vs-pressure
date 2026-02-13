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

analysis_depth = float(sys.argv[1])     
analysis_depth_dz = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

# model names
name_weak_bothfree  = "2D_compositional_subd_lower-res_new_50plates"
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

# text files
text_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

# load files
ref_bothfree   = np.loadtxt((text_ref_bothfree)) 
ref_fixedSP    = np.loadtxt((text_ref_fixedSP))
ref_fixedOP    = np.loadtxt((text_ref_fixedOP))
ref2_bothfree   = np.loadtxt((text_ref2_bothfree)) 
ref2_fixedSP    = np.loadtxt((text_ref2_fixedSP))
ref2_fixedOP    = np.loadtxt((text_ref2_fixedOP))
weak_bothfree  = np.loadtxt((text_weak_bothfree)) 
weak_fixedSP   = np.loadtxt((text_weak_fixedSP))
weak_fixedOP   = np.loadtxt((text_weak_fixedOP))
strong_bothfree = np.loadtxt((text_strong_bothfree)) 
strong_fixedSP     = np.loadtxt((text_strong_fixedSP))
strong_fixedOP     = np.loadtxt((text_strong_fixedOP))

plot_name_png = ''.join(['plots/DP-comparisons/compilations/Visc-vs-K.w-shear-correction.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/Visc-vs-K.w-shear-correction.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])

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


######################################################################

ax=fig.add_subplot(gs[0,0])

# curvature vs. (extracted - analytical - shear stress)
plt.scatter(250*np.ones((len(ref_bothfree),1)), ref_bothfree[:,12]*1000,   color='peru',s=10,edgecolor='black',lw=0.15,zorder=3) 
plt.scatter(250*np.ones((len(ref_fixedSP),1)), ref_fixedSP[:,12]*1000,    color='peru',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
plt.scatter(250*np.ones((len(ref_fixedOP),1)), ref_fixedOP[:,12]*1000,    color='peru',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3)

plt.scatter(500*np.ones((len(ref2_bothfree),1)), ref2_bothfree[:,12]*1000,  color='brown',s=10,edgecolor='black',lw=0.15,zorder=3) 
plt.scatter(500*np.ones((len(ref2_fixedSP),1)), ref2_fixedSP[:,12]*1000,   color='brown',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
plt.scatter(500*np.ones((len(ref2_fixedOP),1)), ref2_fixedOP[:,12]*1000,   color='brown',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3)

plt.scatter(50*np.ones((len(weak_bothfree),1)), weak_bothfree[:,12]*1000,  color='tan',s=10,edgecolor='black',lw=0.15,zorder=3) 
plt.scatter(50*np.ones((len(weak_fixedSP),1)), weak_fixedSP[:,12]*1000,   color='tan',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
plt.scatter(50*np.ones((len(weak_fixedOP),1)), weak_fixedOP[:,12]*1000,   color='tan',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 

plt.scatter(1000*np.ones((len(strong_bothfree),1)), strong_bothfree[:,12]*1000, color='black',s=10,edgecolor='black',lw=0.15,zorder=3) 
plt.scatter(1000*np.ones((len(strong_fixedSP),1)), strong_fixedSP[:,12]*1000,  color='black',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
plt.scatter(1000*np.ones((len(strong_fixedOP),1)), strong_fixedOP[:,12]*1000,  color='black',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 

# plt.ylim(-25, 50); 
# plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("curvature [1/km]",size=6)
plt.ylabel("misfit  [MPa]",size=6)
plt.ylim(-0.001,  0.005)
plt.xlim(0,  1100)

fixed_aspect_ratio(1)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


