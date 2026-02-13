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

name_ref_bothfree 	= "2D_compositional_subd_lower-res_new"
name_ref_fixedSP  	= "2D_compositional_subd_FixedSP_lower-res_new"
name_ref_fixedOP  	= "2D_compositional_subd_FixedOP_lower-res_new"
name_weak_bothfree 	="2D_compositional_subd_lower-res_new_WeakPlates"
name_weak_fixedSP    ="2D_compositional_subd_lower-res_new_FixedSP_WeakPlates"
name_weak_fixedOP    ="2D_compositional_subd_lower-res_new_FixedOP_WeakPlates"
name_strong_bothfree ="2D_compositional_subd_lower-res_new_StiffPlates"
name_strong_fixedSP  ="2D_compositional_subd_lower-res_new_FixedSP_StiffPlates"
name_strong_fixedOP  ="2D_compositional_subd_lower-res_new_FixedOP_StiffPlates"

text_ref_bothfree 	= ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref_fixedSP  	= ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_ref_fixedOP  	= ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_bothfree	= ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_fixedSP  	= ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_weak_fixedOP  	= ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_strong_fixedOP	= ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

ref_bothfree 	= np.loadtxt((text_ref_bothfree)) 
ref_fixedSP  	= np.loadtxt((text_ref_fixedSP))
ref_fixedOP  	= np.loadtxt((text_ref_fixedOP))
weak_bothfree 	= np.loadtxt((text_weak_bothfree)) 
weak_fixedSP  	= np.loadtxt((text_weak_fixedSP))
weak_fixedOP  	= np.loadtxt((text_weak_fixedOP))
strong_bothfree = np.loadtxt((text_strong_bothfree)) 
strong_fixedSP 	= np.loadtxt((text_strong_fixedSP))
strong_fixedOP 	= np.loadtxt((text_strong_fixedOP))

plot_name_png = ''.join(['plots/DP-comparisons/compilations/misfit-vs-curvature_all-mods_z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-bending.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/misfit-vs-curvature_all-mods_z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.w-bending.pdf'])

fig=plt.figure()
gs=GridSpec(3,1) 

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


# misfit vs. curvature
ax=fig.add_subplot(gs[0,0])
# ref, free
plt.scatter(np.abs(ref_bothfree[:,12]),((ref_bothfree[:,3]-ref_bothfree[:,6])-ref_bothfree[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=4) 
# ref, fixed SP
plt.scatter(np.abs(ref_fixedSP[:,12]),((ref_fixedSP[:,3]-ref_fixedSP[:,6])-ref_fixedSP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=4) 
# ref, fixed OP
plt.scatter(np.abs(ref_fixedOP[:,12]),((ref_fixedOP[:,3]-ref_fixedOP[:,6])-ref_fixedOP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=4) 

# weak, free
plt.scatter(np.abs(weak_bothfree[:,12]),((weak_bothfree[:,3]-weak_bothfree[:,6])-weak_bothfree[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=4) 
# weak, fixed SP
plt.scatter(np.abs(weak_fixedSP[:,12]),((weak_fixedSP[:,3]-weak_fixedSP[:,6])-weak_fixedSP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=4) 
# weak, fixed OP
plt.scatter(np.abs(weak_fixedOP[:,12]),((weak_fixedOP[:,3]-weak_fixedOP[:,6])-weak_fixedOP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=4) 

# strong, free
plt.scatter(np.abs(strong_bothfree[:,12]),((strong_bothfree[:,3]-strong_bothfree[:,6])-strong_bothfree[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=4) 
# strong, fixed SP
plt.scatter(np.abs(strong_fixedSP[:,12]),((strong_fixedSP[:,3]-strong_fixedSP[:,6])-strong_fixedSP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=4) 
# strong, fixed OP
plt.scatter(np.abs(strong_fixedOP[:,12]),((strong_fixedOP[:,3]-strong_fixedOP[:,6])-strong_fixedOP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=4) 


# axis stuff
# plt.xlim(0, 300); 
plt.ylim(-25, 50); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("K [1/m]",size=6)
plt.ylabel(r"extracted-analytical" "\n" r" corrected  [MPa]",size=6)


ax=fig.add_subplot(gs[1,0])
# ref, free
plt.scatter(np.abs(ref_bothfree[:,12]),(ref_bothfree[:,3]-ref_bothfree[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed SP
plt.scatter(np.abs(ref_fixedSP[:,12]),(ref_fixedSP[:,3]-ref_fixedSP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed OP
plt.scatter(np.abs(ref_fixedOP[:,12]),(ref_fixedOP[:,3]-ref_fixedOP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 

# weak, free
plt.scatter(np.abs(weak_bothfree[:,12]),(weak_bothfree[:,3]-weak_bothfree[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed SP
plt.scatter(np.abs(weak_fixedSP[:,12]),(weak_fixedSP[:,3]-weak_fixedSP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed OP
plt.scatter(np.abs(weak_fixedOP[:,12]),(weak_fixedOP[:,3]-weak_fixedOP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 

# strong, free
plt.scatter(np.abs(strong_bothfree[:,12]),(strong_bothfree[:,3]-strong_bothfree[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed SP
plt.scatter(np.abs(strong_fixedSP[:,12]),(strong_fixedSP[:,3]-strong_fixedSP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed OP
plt.scatter(np.abs(strong_fixedOP[:,12]),(strong_fixedOP[:,3]-strong_fixedOP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# axis stuff
plt.ylim(-25, 50); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("K [1/m]",size=6)
plt.ylabel(r"extracted-analytical" "\n" r" uncorrected  [MPa]",size=6)


ax=fig.add_subplot(gs[2,0])

# ref, free
plt.scatter(np.abs(ref_bothfree[:,12]),ref_bothfree[:,6]/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=4) 
# ref, fixed SP
plt.scatter(np.abs(ref_fixedSP[:,12]),ref_fixedSP[:,6]/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=4) 
# ref, fixed OP
plt.scatter(np.abs(ref_fixedOP[:,12]),ref_fixedOP[:,6]/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=4) 

# weak, free
plt.scatter(np.abs(weak_bothfree[:,12]),weak_bothfree[:,6]/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=4) 
# weak, fixed SP
plt.scatter(np.abs(weak_fixedSP[:,12]),weak_fixedSP[:,6]/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=4) 
# weak, fixed OP
plt.scatter(np.abs(weak_fixedOP[:,12]),weak_fixedOP[:,6]/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=4) 

# strong, free
plt.scatter(np.abs(strong_bothfree[:,12]),strong_bothfree[:,6]/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=4) 
# strong, fixed SP
plt.scatter(np.abs(strong_fixedSP[:,12]),strong_fixedSP[:,6]/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=4) 
# strong, fixed OP
plt.scatter(np.abs(strong_fixedOP[:,12]),strong_fixedOP[:,6]/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=4) 

# axis stuff
plt.ylim(-50, 25); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("K [1/m]",size=6)
plt.ylabel("shear term  [MPa]",size=6)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
# plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


