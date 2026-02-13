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

analysis_depth = float(sys.argv[1])     # m (depth for DP extraction and central point of shear stress derivative)
analysis_depth_dz = float(sys.argv[2])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

name_ref_bothfree   = "2D_compositional_subd_lower-res_new"
name_ref_fixedSP    = "2D_compositional_subd_FixedSP_lower-res_new"
name_ref_fixedOP    = "2D_compositional_subd_FixedOP_lower-res_new"
name_weak_bothfree  ="2D_compositional_subd_lower-res_new_WeakPlates"
name_weak_fixedSP    ="2D_compositional_subd_lower-res_new_FixedSP_WeakPlates"
name_weak_fixedOP    ="2D_compositional_subd_lower-res_new_FixedOP_WeakPlates"
name_strong_bothfree ="2D_compositional_subd_lower-res_new_StiffPlates"
name_strong_fixedSP  ="2D_compositional_subd_lower-res_new_FixedSP_StiffPlates"
name_strong_fixedOP  ="2D_compositional_subd_lower-res_new_FixedOP_StiffPlates"

text_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

ref_bothfree   = np.loadtxt((text_ref_bothfree)) 
ref_fixedSP    = np.loadtxt((text_ref_fixedSP))
ref_fixedOP    = np.loadtxt((text_ref_fixedOP))
weak_bothfree  = np.loadtxt((text_weak_bothfree)) 
weak_fixedSP   = np.loadtxt((text_weak_fixedSP))
weak_fixedOP   = np.loadtxt((text_weak_fixedOP))
strong_bothfree = np.loadtxt((text_strong_bothfree)) 
strong_fixedSP     = np.loadtxt((text_strong_fixedSP))
strong_fixedOP     = np.loadtxt((text_strong_fixedOP))


plot_name_png = ''.join(['plots/DP-comparisons/compilations/all-mods_z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.wo-shear-correction.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/all-mods_z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.wo-shear-correction.pdf'])

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
curvature_thresh = 1e9
# free plates
for i in range(3,len(weak_bothfree)):
    if np.abs(weak_bothfree[i,12])*1000 < curvature_thresh: 
        plt.scatter(weak_bothfree[i,4]/1.e6,weak_bothfree[i,3]/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(weak_bothfree[i,4]/1.e6,(weak_bothfree[i,3]-weak_bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
# fixed SP
for i in range(3,len(weak_fixedSP)):
    if np.abs(weak_fixedSP[i,12])*1000 < curvature_thresh:     
        plt.scatter(weak_fixedSP[i,4]/1.e6,weak_fixedSP[i,3]/1.e6,color='slateblue',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(weak_fixedSP[i,4]/1.e6,(weak_fixedSP[i,3]-weak_fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
# fixed OP
for i in range(3,len(weak_fixedOP)):
    if np.abs(weak_fixedOP[i,12])*1000 < curvature_thresh: 
        plt.scatter(weak_fixedOP[i,4]/1.e6,weak_fixedOP[i,3]/1.e6,color='darkseagreen',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(weak_fixedOP[i,4]/1.e6,(weak_fixedOP[i,3]-weak_fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected

# axis stuff
plt.xlim(-5,  35); plt.ylim(-5,  35)
plt.plot([-5, 35], [-5, 35], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("extracted stress  [MPa]",size=7)
plt.xlabel("slab buoyancy  [MPa]",size=7)
ax.annotate('$\eta$` = 100', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7.5,color='k') 
fixed_aspect_ratio(1)

ax=fig.add_subplot(gs[0,1])
# free plates
for i in range(3,len(ref_bothfree)):
    if np.abs(ref_bothfree[i,12])*1000 < curvature_thresh: 
        plt.scatter(ref_bothfree[i,4]/1.e6,ref_bothfree[i,3]/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(ref_bothfree[i,4]/1.e6,(ref_bothfree[i,3]-ref_bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
# fixed SP
for i in range(3,len(ref_fixedSP)):
    if np.abs(ref_fixedSP[i,12])*1000 < curvature_thresh:     
        plt.scatter(ref_fixedSP[i,4]/1.e6,ref_fixedSP[i,3]/1.e6,color='slateblue',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(ref_fixedSP[i,4]/1.e6,(ref_fixedSP[i,3]-ref_fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
# fixed OP
for i in range(3,len(ref_fixedOP)):
    if np.abs(ref_fixedOP[i,12])*1000 < curvature_thresh: 
        plt.scatter(ref_fixedOP[i,4]/1.e6,ref_fixedOP[i,3]/1.e6,color='darkseagreen',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(ref_fixedOP[i,4]/1.e6,(ref_fixedOP[i,3]-ref_fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected

# axis stuff
plt.xlim(-5,  35); plt.ylim(-5,  35)
plt.plot([-5, 35], [-5, 35], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab buoyancy  [MPa]",size=7)
ax.annotate('$\eta$` = 500', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7.5,color='k') 
fixed_aspect_ratio(1)

ax=fig.add_subplot(gs[0,2])
# free plates
for i in range(3,len(strong_bothfree)):
    if np.abs(strong_bothfree[i,12])*1000 < curvature_thresh: 
        plt.scatter(strong_bothfree[i,4]/1.e6,strong_bothfree[i,3]/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(strong_bothfree[i,4]/1.e6,(strong_bothfree[i,3]-strong_bothfree[i,6])/1.e6,color='slategray',s=10,edgecolor='black',label='both free',lw=0.25,zorder=4) # corrected
# fixed SP
for i in range(3,len(strong_fixedSP)):
    if np.abs(strong_fixedSP[i,12])*1000 < curvature_thresh:     
        plt.scatter(strong_fixedSP[i,4]/1.e6,strong_fixedSP[i,3]/1.e6,color='slateblue',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(strong_fixedSP[i,4]/1.e6,(strong_fixedSP[i,3]-strong_fixedSP[i,6])/1.e6,color='slateblue',s=10,edgecolor='black',label='fixed SP',lw=0.25,zorder=4) # corrected
# fixed OP
for i in range(3,len(strong_fixedOP)):
    if np.abs(strong_fixedOP[i,12])*1000 < curvature_thresh: 
        plt.scatter(strong_fixedOP[i,4]/1.e6,strong_fixedOP[i,3]/1.e6,color='darkseagreen',s=10,edgecolor='black',lw=0.25,zorder=3) # uncorrected
        # plt.scatter(strong_fixedOP[i,4]/1.e6,(strong_fixedOP[i,3]-strong_fixedOP[i,6])/1.e6,color='darkseagreen',s=10,edgecolor='black',label='fixed OP',lw=0.25,zorder=4) # corrected

# axis stuff
plt.xlim(-5,  35); plt.ylim(-5,  35)
plt.plot([-5, 35], [-5, 35], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab buoyancy  [MPa]",size=7)
ax.annotate('$\eta$` = 2000', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7.5,color='k') 
fixed_aspect_ratio(1)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


