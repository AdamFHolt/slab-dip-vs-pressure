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

analysis_depth_dz = float(sys.argv[1])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[2])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[3])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

analysis_depth1 = 230e3
analysis_depth2 = 330e3
analysis_depth3 = 430e3

name_ref_bothfree 	= "2D_compositional_subd_lower-res_new"
name_ref_fixedSP  	= "2D_compositional_subd_FixedSP_lower-res_new"
name_ref_fixedOP  	= "2D_compositional_subd_FixedOP_lower-res_new"
name_weak_bothfree 	="2D_compositional_subd_lower-res_new_WeakPlates"
name_weak_fixedSP    ="2D_compositional_subd_lower-res_new_FixedSP_WeakPlates"
name_weak_fixedOP    ="2D_compositional_subd_lower-res_new_FixedOP_WeakPlates"
name_strong_bothfree ="2D_compositional_subd_lower-res_new_StiffPlates"
name_strong_fixedSP  ="2D_compositional_subd_lower-res_new_FixedSP_StiffPlates"
name_strong_fixedOP  ="2D_compositional_subd_lower-res_new_FixedOP_StiffPlates"

# shall depths
text1_ref_bothfree 	= ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedSP  	= ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_ref_fixedOP  	= ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_bothfree	= ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedSP  	= ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_weak_fixedOP  	= ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_strong_fixedOP	= ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# intermediate 
text2_ref_bothfree 	= ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedSP  	= ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_ref_fixedOP  	= ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_bothfree	= ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedSP  	= ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_weak_fixedOP  	= ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_strong_fixedOP	= ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# deep
text3_ref_bothfree 	= ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedSP  	= ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_ref_fixedOP  	= ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_bothfree	= ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedSP  	= ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_weak_fixedOP  	= ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_strong_fixedOP	= ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


ref1_bothfree 	= np.loadtxt((text1_ref_bothfree)) 
ref1_fixedSP  	= np.loadtxt((text1_ref_fixedSP))
ref1_fixedOP  	= np.loadtxt((text1_ref_fixedOP))
weak1_bothfree 	= np.loadtxt((text1_weak_bothfree)) 
weak1_fixedSP  	= np.loadtxt((text1_weak_fixedSP))
weak1_fixedOP  	= np.loadtxt((text1_weak_fixedOP))
strong1_bothfree = np.loadtxt((text1_strong_bothfree)) 
strong1_fixedSP 	= np.loadtxt((text1_strong_fixedSP))
strong1_fixedOP 	= np.loadtxt((text1_strong_fixedOP))

ref2_bothfree 	= np.loadtxt((text2_ref_bothfree)) 
ref2_fixedSP  	= np.loadtxt((text2_ref_fixedSP))
ref2_fixedOP  	= np.loadtxt((text2_ref_fixedOP))
weak2_bothfree 	= np.loadtxt((text2_weak_bothfree)) 
weak2_fixedSP  	= np.loadtxt((text2_weak_fixedSP))
weak2_fixedOP  	= np.loadtxt((text2_weak_fixedOP))
strong2_bothfree = np.loadtxt((text2_strong_bothfree)) 
strong2_fixedSP 	= np.loadtxt((text2_strong_fixedSP))
strong2_fixedOP 	= np.loadtxt((text2_strong_fixedOP))

ref3_bothfree 	= np.loadtxt((text3_ref_bothfree)) 
ref3_fixedSP  	= np.loadtxt((text3_ref_fixedSP))
ref3_fixedOP  	= np.loadtxt((text3_ref_fixedOP))
weak3_bothfree 	= np.loadtxt((text3_weak_bothfree)) 
weak3_fixedSP  	= np.loadtxt((text3_weak_fixedSP))
weak3_fixedOP  	= np.loadtxt((text3_weak_fixedOP))
strong3_bothfree = np.loadtxt((text3_strong_bothfree)) 
strong3_fixedSP 	= np.loadtxt((text3_strong_fixedSP))
strong3_fixedOP 	= np.loadtxt((text3_strong_fixedOP))


plot_name_png = ''.join(['plots/DP-comparisons/compilations/misfit-vs-shear-stress-and-curvature_all-mods.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/misfit-vs-shear-stress-and-curvature_all-mods.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])

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

# ref, free
plt.scatter(ref1_bothfree[:,6]/1.e6,(ref1_bothfree[:,3]-ref1_bothfree[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed SP
plt.scatter(ref1_fixedSP[:,6]/1.e6,(ref1_fixedSP[:,3]-ref1_fixedSP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed OP
plt.scatter(ref1_fixedOP[:,6]/1.e6,(ref1_fixedOP[:,3]-ref1_fixedOP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 

# weak, free
plt.scatter(weak1_bothfree[:,6]/1.e6,(weak1_bothfree[:,3]-weak1_bothfree[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed SP
plt.scatter(weak1_fixedSP[:,6]/1.e6,(weak1_fixedSP[:,3]-weak1_fixedSP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed OP
plt.scatter(weak1_fixedOP[:,6]/1.e6,(weak1_fixedOP[:,3]-weak1_fixedOP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 

# strong, free
plt.scatter(strong1_bothfree[:,6]/1.e6,(strong1_bothfree[:,3]-strong1_bothfree[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed SP
plt.scatter(strong1_fixedSP[:,6]/1.e6,(strong1_fixedSP[:,3]-strong1_fixedSP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed OP
plt.scatter(strong1_fixedOP[:,6]/1.e6,(strong1_fixedOP[:,3]-strong1_fixedOP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# axis stuff
# plt.ylim(-25, 50); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("shear stress [MPa]",size=6)
plt.ylabel("extracted-analytical [MPa]",size=6)
plt.xlim(-50,  10); plt.ylim(-30,  30)
plt.plot([-100, 100], [-100, 100], color='black', linewidth=1, zorder=1)
fixed_aspect_ratio(1)
plt.annotate('230 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


# misfit vs. curvature
ax=fig.add_subplot(gs[0,1])

# ref, free
plt.scatter(ref2_bothfree[:,6]/1.e6,(ref2_bothfree[:,3]-ref2_bothfree[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed SP
plt.scatter(ref2_fixedSP[:,6]/1.e6,(ref2_fixedSP[:,3]-ref2_fixedSP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed OP
plt.scatter(ref2_fixedOP[:,6]/1.e6,(ref2_fixedOP[:,3]-ref2_fixedOP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 

# weak, free
plt.scatter(weak2_bothfree[:,6]/1.e6,(weak2_bothfree[:,3]-weak2_bothfree[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed SP
plt.scatter(weak2_fixedSP[:,6]/1.e6,(weak2_fixedSP[:,3]-weak2_fixedSP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed OP
plt.scatter(weak2_fixedOP[:,6]/1.e6,(weak2_fixedOP[:,3]-weak2_fixedOP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 

# strong, free
plt.scatter(strong2_bothfree[:,6]/1.e6,(strong2_bothfree[:,3]-strong2_bothfree[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed SP
plt.scatter(strong2_fixedSP[:,6]/1.e6,(strong2_fixedSP[:,3]-strong2_fixedSP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed OP
plt.scatter(strong2_fixedOP[:,6]/1.e6,(strong2_fixedOP[:,3]-strong2_fixedOP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# axis stuff
# plt.ylim(-25, 50); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("shear stress [MPa]",size=6)
# plt.ylabel("extracted-analytical [MPa]",size=6)
plt.xlim(-40,  20); plt.ylim(-30,  30)
plt.plot([-100, 100], [-100, 100], color='black', linewidth=1, zorder=1)
fixed_aspect_ratio(1)
plt.annotate('330 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


# misfit vs. curvature
ax=fig.add_subplot(gs[0,2])

# ref, free
plt.scatter(ref3_bothfree[:,6]/1.e6,(ref3_bothfree[:,3]-ref3_bothfree[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed SP
plt.scatter(ref3_fixedSP[:,6]/1.e6,(ref3_fixedSP[:,3]-ref3_fixedSP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
# ref, fixed OP
plt.scatter(ref3_fixedOP[:,6]/1.e6,(ref3_fixedOP[:,3]-ref3_fixedOP[:,4])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 

# weak, free
plt.scatter(weak3_bothfree[:,6]/1.e6,(weak3_bothfree[:,3]-weak3_bothfree[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed SP
plt.scatter(weak3_fixedSP[:,6]/1.e6,(weak3_fixedSP[:,3]-weak3_fixedSP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
# weak, fixed OP
plt.scatter(weak3_fixedOP[:,6]/1.e6,(weak3_fixedOP[:,3]-weak3_fixedOP[:,4])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 

# strong, free
plt.scatter(strong3_bothfree[:,6]/1.e6,(strong3_bothfree[:,3]-strong3_bothfree[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed SP
plt.scatter(strong3_fixedSP[:,6]/1.e6,(strong3_fixedSP[:,3]-strong3_fixedSP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
# strong, fixed OP
plt.scatter(strong3_fixedOP[:,6]/1.e6,(strong3_fixedOP[:,3]-strong3_fixedOP[:,4])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# axis stuff
# plt.ylim(-25, 50); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("shear stress [MPa]",size=6)
# plt.ylabel("extracted-analytical [MPa]",size=6)
plt.xlim(-35,  25); plt.ylim(-30,  30)
plt.plot([-100, 100], [-100, 100], color='black', linewidth=1, zorder=1)
fixed_aspect_ratio(1)
plt.annotate('430 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


######################################################################
######################################################################

ax=fig.add_subplot(gs[1,0])
# curvature vs. (extracted - analytical - shear stress)
plt.scatter(ref1_bothfree[:,12]*1000,(ref1_bothfree[:,3]-ref1_bothfree[:,4]-ref1_bothfree[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(ref1_fixedSP[:,12]*1000,(ref1_fixedSP[:,3]-ref1_fixedSP[:,4]-ref1_fixedSP[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(ref1_fixedOP[:,12]*1000,(ref1_fixedOP[:,3]-ref1_fixedOP[:,4]-ref1_fixedOP[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak1_bothfree[:,12]*1000,(weak1_bothfree[:,3]-weak1_bothfree[:,4]-weak1_bothfree[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak1_fixedSP[:,12]*1000,(weak1_fixedSP[:,3]-weak1_fixedSP[:,4]-weak1_fixedSP[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak1_fixedOP[:,12]*1000,(weak1_fixedOP[:,3]-weak1_fixedOP[:,4]-weak1_fixedOP[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong1_bothfree[:,12]*1000,(strong1_bothfree[:,3]-strong1_bothfree[:,4]-strong1_bothfree[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong1_fixedSP[:,12]*1000,(strong1_fixedSP[:,3]-strong1_fixedSP[:,4]-strong1_fixedSP[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong1_fixedOP[:,12]*1000,(strong1_fixedOP[:,3]-strong1_fixedOP[:,4]-strong1_fixedOP[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# plt.ylim(-25, 50); 
# plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.xlabel("curvature [1/km]",size=6)
plt.ylabel("extracted-analytical-shear [MPa]",size=6)
# plt.xlim(-50,  10)
plt.ylim(-20,  55)
# plt.plot([-100, 100], [-100, 100], color='black', linewidth=1, zorder=1)
fixed_aspect_ratio(1)
plt.annotate('230 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')

ax=fig.add_subplot(gs[1,1])
# curvature vs. (extracted - analytical - shear stress)
plt.scatter(ref2_bothfree[:,12]*1000,(ref2_bothfree[:,3]-ref2_bothfree[:,4]-ref2_bothfree[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(ref2_fixedSP[:,12]*1000,(ref2_fixedSP[:,3]-ref2_fixedSP[:,4]-ref2_fixedSP[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(ref2_fixedOP[:,12]*1000,(ref2_fixedOP[:,3]-ref2_fixedOP[:,4]-ref2_fixedOP[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak2_bothfree[:,12]*1000,(weak2_bothfree[:,3]-weak2_bothfree[:,4]-weak2_bothfree[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak2_fixedSP[:,12]*1000,(weak2_fixedSP[:,3]-weak2_fixedSP[:,4]-weak2_fixedSP[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak2_fixedOP[:,12]*1000,(weak2_fixedOP[:,3]-weak2_fixedOP[:,4]-weak2_fixedOP[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong2_bothfree[:,12]*1000,(strong2_bothfree[:,3]-strong2_bothfree[:,4]-strong2_bothfree[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong2_fixedSP[:,12]*1000,(strong2_fixedSP[:,3]-strong2_fixedSP[:,4]-strong2_fixedSP[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong2_fixedOP[:,12]*1000,(strong2_fixedOP[:,3]-strong2_fixedOP[:,4]-strong2_fixedOP[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# plt.ylim(-25, 50); 
# plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("curvature [1/km]",size=6)
# plt.xlim(-50,  10); 
plt.ylim(-20,  55)
# plt.plot([-100, 100], [-100, 100], color='black', linewidth=1, zorder=1)
fixed_aspect_ratio(1)
plt.annotate('330 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')

ax=fig.add_subplot(gs[1,2])
# curvature vs. (extracted - analytical - shear stress)
plt.scatter(ref3_bothfree[:,12]*1000,(ref3_bothfree[:,3]-ref3_bothfree[:,4]-ref3_bothfree[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(ref3_fixedSP[:,12]*1000,(ref3_fixedSP[:,3]-ref3_fixedSP[:,4]-ref3_fixedSP[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(ref3_fixedOP[:,12]*1000,(ref3_fixedOP[:,3]-ref3_fixedOP[:,4]-ref3_fixedOP[:,6])/1.e6,color='slategray',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak3_bothfree[:,12]*1000,(weak3_bothfree[:,3]-weak3_bothfree[:,4]-weak3_bothfree[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak3_fixedSP[:,12]*1000,(weak3_fixedSP[:,3]-weak3_fixedSP[:,4]-weak3_fixedSP[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(weak3_fixedOP[:,12]*1000,(weak3_fixedOP[:,3]-weak3_fixedOP[:,4]-weak3_fixedOP[:,6])/1.e6,color='blue',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong3_bothfree[:,12]*1000,(strong3_bothfree[:,3]-strong3_bothfree[:,4]-strong3_bothfree[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong3_fixedSP[:,12]*1000,(strong3_fixedSP[:,3]-strong3_fixedSP[:,4]-strong3_fixedSP[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 
plt.scatter(strong3_fixedOP[:,12]*1000,(strong3_fixedOP[:,3]-strong3_fixedOP[:,4]-strong3_fixedOP[:,6])/1.e6,color='red',s=10,edgecolor='black',lw=0.25,zorder=3) 

# plt.ylim(-25, 50); 
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("curvature [1/km]",size=6)
# plt.xlim(-50,  10); 
plt.ylim(-20,  55)
# plt.plot([-100, 100], [-100, 100], color='black', linewidth=1, zorder=1)
fixed_aspect_ratio(1)
plt.annotate('430 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


