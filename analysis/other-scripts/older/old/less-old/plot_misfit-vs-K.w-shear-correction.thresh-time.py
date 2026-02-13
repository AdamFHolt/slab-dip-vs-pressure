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
analysis_depth_dz = float(sys.argv[2])  # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)
threshold_time = 30.                    # Myrs

# model names
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
text_new_bothfree=''.join(['text_files/new/',name_new_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_new_fixedSP= ''.join(['text_files/new/',name_new_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text_new_fixedOP= ''.join(['text_files/new/',name_new_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

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
new_bothfree = np.loadtxt((text_new_bothfree)) 
new_fixedSP     = np.loadtxt((text_new_fixedSP))
new_fixedOP     = np.loadtxt((text_new_fixedOP))

plot_name_png = ''.join(['plots/DP-comparisons/compilations/misfit-vs-K.w-shear-correction.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.thresh-time',str(threshold_time),'Myr.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/misfit-vs-K.w-shear-correction.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.thresh-time',str(threshold_time),'Myr.pdf'])

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

for i in range(len(weak_bothfree)):
	if weak_bothfree[i,0] > threshold_time:
		plt.scatter(weak_bothfree[i,4]/1.e6,  (weak_bothfree[i,3]-weak_bothfree[i,6])/1.e6,    color='tan',     s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(weak_bothfree[i,4]/1.e6,  (weak_bothfree[i,3]-weak_bothfree[i,6])/1.e6,    color='tan',     s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(weak_fixedSP)):
	if weak_fixedSP[i,0] > threshold_time:
		plt.scatter(weak_fixedSP[i,4]/1.e6,   (weak_fixedSP[i,3]-weak_fixedSP[i,6])/1.e6,      color='tan',     s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
	else:
		plt.scatter(weak_fixedSP[i,4]/1.e6,   (weak_fixedSP[i,3]-weak_fixedSP[i,6])/1.e6,      color='tan',     s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(weak_fixedOP)):
	if weak_fixedOP[i,0] > threshold_time:
		plt.scatter(weak_fixedOP[i,4]/1.e6,   (weak_fixedOP[i,3]-weak_fixedOP[i,6])/1.e6,      color='tan',     s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(weak_fixedOP[i,4]/1.e6,   (weak_fixedOP[i,3]-weak_fixedOP[i,6])/1.e6,      color='tan',     s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(ref_bothfree)):
	if ref_bothfree[i,0] > threshold_time:
		plt.scatter(ref_bothfree[i,4]/1.e6,   (ref_bothfree[i,3]-ref_bothfree[i,6])/1.e6,      color='peru',s=10,edgecolor='black',lw=0.15,zorder=3)
	else:
		plt.scatter(ref_bothfree[i,4]/1.e6,   (ref_bothfree[i,3]-ref_bothfree[i,6])/1.e6,      color='peru',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2)
for i in range(len(ref_fixedSP)):
	if ref_fixedSP[i,0] > threshold_time: 
		plt.scatter(ref_fixedSP[i,4]/1.e6,    (ref_fixedSP[i,3]-ref_fixedSP[i,6])/1.e6,   	   color='peru',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3)
	else:
		plt.scatter(ref_fixedSP[i,4]/1.e6,    (ref_fixedSP[i,3]-ref_fixedSP[i,6])/1.e6,   	   color='peru',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2)
for i in range(len(ref_fixedOP)):
	if ref_fixedOP[i,0] > threshold_time: 
		plt.scatter(ref_fixedOP[i,4]/1.e6,    (ref_fixedOP[i,3]-ref_fixedOP[i,6])/1.e6,        color='peru',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref_fixedOP[i,4]/1.e6,    (ref_fixedOP[i,3]-ref_fixedOP[i,6])/1.e6,        color='peru',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(new_bothfree)):
	if new_bothfree[i,0] > threshold_time:
		plt.scatter(new_bothfree[i,4]/1.e6,   (new_bothfree[i,3]-new_bothfree[i,6])/1.e6,      color='firebrick',s=10,edgecolor='black',lw=0.15,zorder=3)
	else: 
		plt.scatter(new_bothfree[i,4]/1.e6,   (new_bothfree[i,3]-new_bothfree[i,6])/1.e6,      color='firebrick',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2)
for i in range(len(new_fixedSP)):
	if new_fixedSP[i,0] > threshold_time: 
		plt.scatter(new_fixedSP[i,4]/1.e6,    (new_fixedSP[i,3]-new_fixedSP[i,6])/1.e6,   	   color='firebrick',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3)
	else:
		plt.scatter(new_fixedSP[i,4]/1.e6,    (new_fixedSP[i,3]-new_fixedSP[i,6])/1.e6,   	   color='firebrick',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2)
for i in range(len(new_fixedOP)):
	if new_fixedOP[i,0] > threshold_time: 
		plt.scatter(new_fixedOP[i,4]/1.e6,    (new_fixedOP[i,3]-new_fixedOP[i,6])/1.e6,        color='firebrick',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(new_fixedOP[i,4]/1.e6,    (new_fixedOP[i,3]-new_fixedOP[i,6])/1.e6,        color='firebrick',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(ref2_bothfree)):
	if ref2_bothfree[i,0] > threshold_time:
		plt.scatter(ref2_bothfree[i,4]/1.e6,   (ref2_bothfree[i,3]-ref2_bothfree[i,6])/1.e6,      color='brown',s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref2_bothfree[i,4]/1.e6,   (ref2_bothfree[i,3]-ref2_bothfree[i,6])/1.e6,      color='brown',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(ref2_fixedSP)):
	if ref2_fixedSP[i,0] > threshold_time:
		plt.scatter(ref2_fixedSP[i,4]/1.e6,    (ref2_fixedSP[i,3]-ref2_fixedSP[i,6])/1.e6,   	   color='brown',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref2_fixedSP[i,4]/1.e6,    (ref2_fixedSP[i,3]-ref2_fixedSP[i,6])/1.e6,   	   color='brown',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(ref2_fixedOP)):
	if ref2_fixedOP[i,0] > threshold_time:
		plt.scatter(ref2_fixedOP[i,4]/1.e6,    (ref2_fixedOP[i,3]-ref2_fixedOP[i,6])/1.e6,        color='brown',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref2_fixedOP[i,4]/1.e6,    (ref2_fixedOP[i,3]-ref2_fixedOP[i,6])/1.e6,        color='brown',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(strong_bothfree)):
	if strong_bothfree[i,0] > threshold_time:
		plt.scatter(strong_bothfree[i,4]/1.e6,(strong_bothfree[i,3]-strong_bothfree[i,6])/1.e6,color='black',      s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(strong_bothfree[i,4]/1.e6,(strong_bothfree[i,3]-strong_bothfree[i,6])/1.e6,color='black',      s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(strong_fixedSP)):
	if strong_fixedSP[i,0] > threshold_time:
		plt.scatter(strong_fixedSP[i,4]/1.e6, (strong_fixedSP[i,3]-strong_fixedSP[i,6])/1.e6,  color='black',      s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
	else:
		plt.scatter(strong_fixedSP[i,4]/1.e6, (strong_fixedSP[i,3]-strong_fixedSP[i,6])/1.e6,  color='black',      s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(strong_fixedOP)):
	if strong_fixedOP[i,0] > threshold_time:		
		plt.scatter(strong_fixedOP[i,4]/1.e6, (strong_fixedOP[i,3]-strong_fixedOP[i,6])/1.e6,  color='black',      s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(strong_fixedOP[i,4]/1.e6, (strong_fixedOP[i,3]-strong_fixedOP[i,6])/1.e6,  color='black',      s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

# axis stuff
plt.xlim(-15,  45); plt.ylim(-15,  45)
plt.plot([-15, 45], [-15, 45], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel(r'$\mathregular{\Delta P + h\partial _s \tau}$  [MPa]',size=7)
plt.xlabel(r'$\mathregular{\Delta \rho g_n h}$  [MPa]',size=7)
ax.set_xticks( [-10,0,10,20,30,40] )
fixed_aspect_ratio(1)


ax=fig.add_subplot(gs[0,1])

# curvature vs. (extracted - analytical - shear stress)
for i in range(len(weak_bothfree)):
	if weak_bothfree[i,0] > threshold_time:
		plt.scatter(weak_bothfree[i,12]*1000,  (weak_bothfree[i,3]-weak_bothfree[i,4]-weak_bothfree[i,6])/1.e6,      color='tan',s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(weak_bothfree[i,12]*1000,  (weak_bothfree[i,3]-weak_bothfree[i,4]-weak_bothfree[i,6])/1.e6,      color='tan',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(weak_fixedSP)):
	if weak_fixedSP[i,0] > threshold_time:
		plt.scatter(weak_fixedSP[i,12]*1000,   (weak_fixedSP[i,3]-weak_fixedSP[i,4]-weak_fixedSP[i,6])/1.e6,         color='tan',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
	else:
		plt.scatter(weak_fixedSP[i,12]*1000,   (weak_fixedSP[i,3]-weak_fixedSP[i,4]-weak_fixedSP[i,6])/1.e6,         color='tan',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(weak_fixedOP)):
	if weak_fixedOP[i,0] > threshold_time:
		plt.scatter(weak_fixedOP[i,12]*1000,   (weak_fixedOP[i,3]-weak_fixedOP[i,4]-weak_fixedOP[i,6])/1.e6,         color='tan',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(weak_fixedOP[i,12]*1000,   (weak_fixedOP[i,3]-weak_fixedOP[i,4]-weak_fixedOP[i,6])/1.e6,         color='tan',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

for i in range(len(ref_bothfree)):
	if ref_bothfree[i,0] > threshold_time:
		plt.scatter(ref_bothfree[i,12]*1000,   (ref_bothfree[i,3]-ref_bothfree[i,4]-ref_bothfree[i,6])/1.e6,         color='peru',s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref_bothfree[i,12]*1000,   (ref_bothfree[i,3]-ref_bothfree[i,4]-ref_bothfree[i,6])/1.e6,         color='peru',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(ref_fixedSP)):
	if ref_fixedSP[i,0] > threshold_time:
		plt.scatter(ref_fixedSP[i,12]*1000,    (ref_fixedSP[i,3]-ref_fixedSP[i,4]-ref_fixedSP[i,6])/1.e6,            color='peru',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3)
	else:
		plt.scatter(ref_fixedSP[i,12]*1000,    (ref_fixedSP[i,3]-ref_fixedSP[i,4]-ref_fixedSP[i,6])/1.e6,            color='peru',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2)
for i in range(len(ref_fixedOP)):
	if ref_fixedOP[i,0] > threshold_time: 
		plt.scatter(ref_fixedOP[i,12]*1000,    (ref_fixedOP[i,3]-ref_fixedOP[i,4]-ref_fixedOP[i,6])/1.e6,            color='peru',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3)
	else:
		plt.scatter(ref_fixedOP[i,12]*1000,    (ref_fixedOP[i,3]-ref_fixedOP[i,4]-ref_fixedOP[i,6])/1.e6,            color='peru',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2)

for i in range(len(new_bothfree)):
	if new_bothfree[i,0] > threshold_time:
		plt.scatter(new_bothfree[i,12]*1000,   (new_bothfree[i,3]-new_bothfree[i,4]-new_bothfree[i,6])/1.e6,         color='firebrick',s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(new_bothfree[i,12]*1000,   (new_bothfree[i,3]-new_bothfree[i,4]-new_bothfree[i,6])/1.e6,         color='firebrick',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(new_fixedSP)):
	if new_fixedSP[i,0] > threshold_time:
		plt.scatter(new_fixedSP[i,12]*1000,    (new_fixedSP[i,3]-new_fixedSP[i,4]-new_fixedSP[i,6])/1.e6,            color='firebrick',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3)
	else:
		plt.scatter(new_fixedSP[i,12]*1000,    (new_fixedSP[i,3]-new_fixedSP[i,4]-new_fixedSP[i,6])/1.e6,            color='firebrick',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2)
for i in range(len(new_fixedOP)):
	if new_fixedOP[i,0] > threshold_time: 
		plt.scatter(new_fixedOP[i,12]*1000,    (new_fixedOP[i,3]-new_fixedOP[i,4]-new_fixedOP[i,6])/1.e6,            color='firebrick',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3)
	else:
		plt.scatter(new_fixedOP[i,12]*1000,    (new_fixedOP[i,3]-new_fixedOP[i,4]-new_fixedOP[i,6])/1.e6,            color='firebrick',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2)

for i in range(len(ref2_bothfree)):
	if ref2_bothfree[i,0] > threshold_time:
		plt.scatter(ref2_bothfree[i,12]*1000,   (ref2_bothfree[i,3]-ref2_bothfree[i,4]-ref2_bothfree[i,6])/1.e6,         color='maroon',s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref2_bothfree[i,12]*1000,   (ref2_bothfree[i,3]-ref2_bothfree[i,4]-ref2_bothfree[i,6])/1.e6,         color='maroon',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(ref2_fixedSP)):
	if ref2_fixedSP[i,0] > threshold_time:
		plt.scatter(ref2_fixedSP[i,12]*1000,    (ref2_fixedSP[i,3]-ref2_fixedSP[i,4]-ref2_fixedSP[i,6])/1.e6,            color='maroon',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
	else:
		plt.scatter(ref2_fixedSP[i,12]*1000,    (ref2_fixedSP[i,3]-ref2_fixedSP[i,4]-ref2_fixedSP[i,6])/1.e6,            color='maroon',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(ref2_fixedOP)):
	if ref2_fixedOP[i,0] > threshold_time:
		plt.scatter(ref2_fixedOP[i,12]*1000,    (ref2_fixedOP[i,3]-ref2_fixedOP[i,4]-ref2_fixedOP[i,6])/1.e6,            color='maroon',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3)
	else:
		plt.scatter(ref2_fixedOP[i,12]*1000,    (ref2_fixedOP[i,3]-ref2_fixedOP[i,4]-ref2_fixedOP[i,6])/1.e6,            color='maroon',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2)

for i in range(len(strong_bothfree)):
	if strong_bothfree[i,0] > threshold_time:
		plt.scatter(strong_bothfree[i,12]*1000,(strong_bothfree[i,3]-strong_bothfree[i,4]-strong_bothfree[i,6])/1.e6,color='black',s=10,edgecolor='black',lw=0.15,zorder=3) 
	else:
		plt.scatter(strong_bothfree[i,12]*1000,(strong_bothfree[i,3]-strong_bothfree[i,4]-strong_bothfree[i,6])/1.e6,color='black',s=5,edgecolor='black',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(strong_fixedSP)):
	if strong_fixedSP[i,0] > threshold_time:
		plt.scatter(strong_fixedSP[i,12]*1000, (strong_fixedSP[i,3]-strong_fixedSP[i,4]-strong_fixedSP[i,6])/1.e6,   color='black',s=10,edgecolor='black',marker='v',lw=0.15,zorder=3) 
	else:
		plt.scatter(strong_fixedSP[i,12]*1000, (strong_fixedSP[i,3]-strong_fixedSP[i,4]-strong_fixedSP[i,6])/1.e6,   color='black',s=5,edgecolor='black',marker='v',lw=0.15,zorder=3,alpha=0.2) 
for i in range(len(strong_fixedOP)):
	if strong_fixedOP[i,0] > threshold_time:
		plt.scatter(strong_fixedOP[i,12]*1000, (strong_fixedOP[i,3]-strong_fixedOP[i,4]-strong_fixedOP[i,6])/1.e6,   color='black',s=10,edgecolor='black',marker='^',lw=0.15,zorder=3) 
	else:
		plt.scatter(strong_fixedOP[i,12]*1000, (strong_fixedOP[i,3]-strong_fixedOP[i,4]-strong_fixedOP[i,6])/1.e6,   color='black',s=5,edgecolor='black',marker='^',lw=0.15,zorder=3,alpha=0.2) 

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("curvature [1/km]",size=6)
plt.ylabel("misfit  [MPa]",size=6)
plt.ylim(-20,  55)
fixed_aspect_ratio(1)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


