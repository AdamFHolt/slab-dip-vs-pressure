#!/bin/python3.10

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
from functions import get_avg_forces_nondim

analysis_depth  = float(sys.argv[1]) 
analysis_depth_dz = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)


plot_name_png = ''.join(['plots/DP-comparisons/compilations/NEW/force-components-averaged-vs-viscosity.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/NEW/force-components-averaged-vs-viscosity.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])

# 50
name1_bothfree 	= "2D_compositional_subd_lower-res_new_50plates"
name1_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name1_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_50plates"
# 250
name3_bothfree 	= "2D_compositional_subd_lower-res_new_250plates"
name3_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_250plates"
name3_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_250plates"
# 500
name4_bothfree 	= "2D_compositional_subd_lower-res_new"
name4_fixedSP  	= "2D_compositional_subd_FixedSP_lower-res_new"
name4_fixedOP  	= "2D_compositional_subd_FixedOP_lower-res_new"
# 1000
name5_bothfree = "2D_compositional_subd_lower-res_new_1000plates"
name5_fixedSP  = "2D_compositional_subd_lower-res_new_FixedSP_1000plates"
name5_fixedOP  = "2D_compositional_subd_lower-res_new_FixedOP_1000plates"
# 375
name7_bothfree 	= "2D_compositional_subd_lower-res_new_375plates"
name7_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_375plates"
name7_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_375plates"

# text files
text50_bothfree 		= ''.join(['text_files/new2/',name1_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedSP  		= ''.join(['text_files/new2/',name1_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedOP  		= ''.join(['text_files/new2/',name1_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_bothfree    	= ''.join(['text_files/new2/',name3_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedSP     	= ''.join(['text_files/new2/',name3_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedOP		    = ''.join(['text_files/new2/',name3_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_bothfree 		= ''.join(['text_files/new2/',name4_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedSP  		= ''.join(['text_files/new2/',name4_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedOP  		= ''.join(['text_files/new2/',name4_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_bothfree		= ''.join(['text_files/new2/',name5_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedSP  		= ''.join(['text_files/new2/',name5_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedOP  		= ''.join(['text_files/new2/',name5_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_bothfree    	= ''.join(['text_files/new2/',name7_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedSP     	= ''.join(['text_files/new2/',name7_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedOP	    	= ''.join(['text_files/new2/',name7_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


# load in models
m50_bothfree 	= np.loadtxt((text50_bothfree)) 
m50_fixedSP  	= np.loadtxt((text50_fixedSP))
m50_fixedOP  	= np.loadtxt((text50_fixedOP))
m250_bothfree 	= np.loadtxt((text250_bothfree)) 
m250_fixedSP 	= np.loadtxt((text250_fixedSP))
m250_fixedOP 	= np.loadtxt((text250_fixedOP))
m500_bothfree 	= np.loadtxt((text500_bothfree)) 
m500_fixedSP  	= np.loadtxt((text500_fixedSP))
m500_fixedOP  	= np.loadtxt((text500_fixedOP))
m1000_bothfree 	= np.loadtxt((text1000_bothfree)) 
m1000_fixedSP  	= np.loadtxt((text1000_fixedSP))
m1000_fixedOP  	= np.loadtxt((text1000_fixedOP))
m375_bothfree 	= np.loadtxt((text375_bothfree)) 
m375_fixedSP 	= np.loadtxt((text375_fixedSP))
m375_fixedOP 	= np.loadtxt((text375_fixedOP))
# --

# col 1 = visc, col 2 = DP anal, col 3 = full mod, col 4 = DP, col 5 = shear, col 6 = norm
fixedSP_array = np.zeros((5,7))
fixedSP_array[0:5,0] = [50,250,375,500,675]
fixedSP_array[0,1], fixedSP_array[0,2], fixedSP_array[0,3], fixedSP_array[0,4], fixedSP_array[0,5], fixedSP_array[0,6] = get_avg_forces_nondim(m50_fixedSP)
fixedSP_array[1,1], fixedSP_array[1,2], fixedSP_array[1,3], fixedSP_array[1,4], fixedSP_array[1,5], fixedSP_array[1,6] = get_avg_forces_nondim(m250_fixedSP)
fixedSP_array[2,1], fixedSP_array[2,2], fixedSP_array[2,3], fixedSP_array[2,4], fixedSP_array[2,5], fixedSP_array[2,6] = get_avg_forces_nondim(m375_fixedSP)
fixedSP_array[3,1], fixedSP_array[3,2], fixedSP_array[3,3], fixedSP_array[3,4], fixedSP_array[3,5], fixedSP_array[3,6] = get_avg_forces_nondim(m500_fixedSP)
fixedSP_array[4,1], fixedSP_array[4,2], fixedSP_array[4,3], fixedSP_array[4,4], fixedSP_array[4,5], fixedSP_array[4,6] = get_avg_forces_nondim(m1000_fixedSP)

bothfree_array = np.zeros((5,7))
bothfree_array[0:5,0] = [50,250,375,500,675]
bothfree_array[0,1], bothfree_array[0,2], bothfree_array[0,3], bothfree_array[0,4], bothfree_array[0,5], bothfree_array[0,6] = get_avg_forces_nondim(m50_bothfree)
bothfree_array[1,1], bothfree_array[1,2], bothfree_array[1,3], bothfree_array[1,4], bothfree_array[1,5], bothfree_array[1,6] = get_avg_forces_nondim(m250_bothfree)
bothfree_array[2,1], bothfree_array[2,2], bothfree_array[2,3], bothfree_array[2,4], bothfree_array[2,5], bothfree_array[2,6] = get_avg_forces_nondim(m375_bothfree)
bothfree_array[3,1], bothfree_array[3,2], bothfree_array[3,3], bothfree_array[3,4], bothfree_array[3,5], bothfree_array[3,6] = get_avg_forces_nondim(m500_bothfree)
bothfree_array[4,1], bothfree_array[4,2], bothfree_array[4,3], bothfree_array[4,4], bothfree_array[4,5], bothfree_array[4,6] = get_avg_forces_nondim(m1000_bothfree)

fixedOP_array = np.zeros((5,7))
fixedOP_array[0:5,0] = [50,250,375,500,675]
fixedOP_array[0,1], fixedOP_array[0,2], fixedOP_array[0,3], fixedOP_array[0,4], fixedOP_array[0,5], fixedOP_array[0,6] = get_avg_forces_nondim(m50_fixedOP)
fixedOP_array[1,1], fixedOP_array[1,2], fixedOP_array[1,3], fixedOP_array[1,4], fixedOP_array[1,5], fixedOP_array[1,6] = get_avg_forces_nondim(m250_fixedOP)
fixedOP_array[2,1], fixedOP_array[2,2], fixedOP_array[2,3], fixedOP_array[2,4], fixedOP_array[2,5], fixedOP_array[2,6] = get_avg_forces_nondim(m375_fixedOP)
fixedOP_array[3,1], fixedOP_array[3,2], fixedOP_array[3,3], fixedOP_array[3,4], fixedOP_array[3,5], fixedOP_array[3,6] = get_avg_forces_nondim(m500_fixedOP)
fixedOP_array[4,1], fixedOP_array[4,2], fixedOP_array[4,3], fixedOP_array[4,4], fixedOP_array[4,5], fixedOP_array[4,6] = get_avg_forces_nondim(m1000_fixedOP)


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


ax=fig.add_subplot(gs[0,0])

# fixed SP models
plt.scatter(fixedSP_array[:,0],fixedSP_array[:,2]/fixedSP_array[:,1],s=20,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')
plt.scatter(fixedSP_array[:,0],fixedSP_array[:,3]/fixedSP_array[:,1],s=20,color='red',edgecolor='red',linewidth=0.25,zorder=3,marker='s')
plt.scatter(fixedSP_array[:,0],fixedSP_array[:,4]/fixedSP_array[:,1],s=20,color='blue',edgecolor='blue',linewidth=0.25,zorder=3,marker='^')
plt.scatter(fixedSP_array[:,0],fixedSP_array[:,5]/fixedSP_array[:,1],s=20,color='green',edgecolor='green',linewidth=0.25,zorder=3,marker='v')


ax.set_xticks( [50,250, 375, 500,675] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.75, zorder=0)
plt.axhline(y=1, color='gray',linestyle='--',linewidth=0.75, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("force/buoyancy",size=6)
plt.ylim(-0.5,  1.5); 

ax=fig.add_subplot(gs[1,0])

# both free models
plt.scatter(bothfree_array[:,0],bothfree_array[:,2]/bothfree_array[:,1],s=20,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')
plt.scatter(bothfree_array[:,0],bothfree_array[:,3]/bothfree_array[:,1],s=20,color='red',edgecolor='red',linewidth=0.25,zorder=3,marker='s')
plt.scatter(bothfree_array[:,0],bothfree_array[:,4]/bothfree_array[:,1],s=20,color='blue',edgecolor='blue',linewidth=0.25,zorder=3,marker='^')
plt.scatter(bothfree_array[:,0],bothfree_array[:,5]/bothfree_array[:,1],s=20,color='green',edgecolor='green',linewidth=0.25,zorder=3,marker='v')

ax.set_xticks( [50,250, 375, 500,675] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.75, zorder=0)
plt.axhline(y=1, color='gray',linestyle='--',linewidth=0.75, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("force/buoyancy",size=6)
plt.ylim(-0.5,  1.5); 


ax=fig.add_subplot(gs[2,0])

# fixed OP models
plt.scatter(fixedOP_array[:,0],fixedOP_array[:,2]/fixedOP_array[:,1],s=20,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')
plt.scatter(fixedOP_array[:,0],fixedOP_array[:,3]/fixedOP_array[:,1],s=20,color='red',edgecolor='red',linewidth=0.25,zorder=3,marker='s')
plt.scatter(fixedOP_array[:,0],fixedOP_array[:,4]/fixedOP_array[:,1],s=20,color='blue',edgecolor='blue',linewidth=0.25,zorder=3,marker='^')
plt.scatter(fixedOP_array[:,0],fixedOP_array[:,5]/fixedOP_array[:,1],s=20,color='green',edgecolor='green',linewidth=0.25,zorder=3,marker='v')


ax.set_xticks( [50,250, 375, 500,675] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.75, zorder=0)
plt.axhline(y=1, color='gray',linestyle='--',linewidth=0.75, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("force/buoyancy",size=6)
plt.ylim(-0.5,  1.5); 


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


