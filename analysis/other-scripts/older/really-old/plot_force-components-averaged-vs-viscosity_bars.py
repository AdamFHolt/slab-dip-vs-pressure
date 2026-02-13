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


plot_name_png = ''.join(['plots/DP-comparisons/compilations/NEW/force-components-averaged-vs-viscosity_bars.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/NEW/force-components-averaged-vs-viscosity_bars.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])

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
text50_bothfree 		= ''.join(['text_files/new2/old/',name1_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedSP  		= ''.join(['text_files/new2/old/',name1_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedOP  		= ''.join(['text_files/new2/old/',name1_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_bothfree    	= ''.join(['text_files/new2/old/',name3_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedSP     	= ''.join(['text_files/new2/old/',name3_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedOP		    = ''.join(['text_files/new2/old/',name3_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_bothfree 		= ''.join(['text_files/new2/old/',name4_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedSP  		= ''.join(['text_files/new2/old/',name4_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedOP  		= ''.join(['text_files/new2/old/',name4_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_bothfree		= ''.join(['text_files/new2/old/',name5_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedSP  		= ''.join(['text_files/new2/old/',name5_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedOP  		= ''.join(['text_files/new2/old/',name5_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_bothfree    	= ''.join(['text_files/new2/old/',name7_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedSP     	= ''.join(['text_files/new2/old/',name7_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedOP	    	= ''.join(['text_files/new2/old/',name7_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


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
fixedSP_array = np.zeros((5,6))
fixedSP_array[0:5,0] = [50,250,375,500,675]
fixedSP_array[0,1], fixedSP_array[0,2], fixedSP_array[0,3], fixedSP_array[0,4], fixedSP_array[0,5] = get_avg_forces_nondim(m50_fixedSP)
fixedSP_array[1,1], fixedSP_array[1,2], fixedSP_array[1,3], fixedSP_array[1,4], fixedSP_array[1,5] = get_avg_forces_nondim(m250_fixedSP)
fixedSP_array[2,1], fixedSP_array[2,2], fixedSP_array[2,3], fixedSP_array[2,4], fixedSP_array[2,5] = get_avg_forces_nondim(m375_fixedSP)
fixedSP_array[3,1], fixedSP_array[3,2], fixedSP_array[3,3], fixedSP_array[3,4], fixedSP_array[3,5] = get_avg_forces_nondim(m500_fixedSP)
fixedSP_array[4,1], fixedSP_array[4,2], fixedSP_array[4,3], fixedSP_array[4,4], fixedSP_array[4,5] = get_avg_forces_nondim(m1000_fixedSP)

bothfree_array = np.zeros((5,6))
bothfree_array[0:5,0] = [50,250,375,500,675]
bothfree_array[0,1], bothfree_array[0,2], bothfree_array[0,3], bothfree_array[0,4], bothfree_array[0,5] = get_avg_forces_nondim(m50_bothfree)
bothfree_array[1,1], bothfree_array[1,2], bothfree_array[1,3], bothfree_array[1,4], bothfree_array[1,5] = get_avg_forces_nondim(m250_bothfree)
bothfree_array[2,1], bothfree_array[2,2], bothfree_array[2,3], bothfree_array[2,4], bothfree_array[2,5] = get_avg_forces_nondim(m375_bothfree)
bothfree_array[3,1], bothfree_array[3,2], bothfree_array[3,3], bothfree_array[3,4], bothfree_array[3,5] = get_avg_forces_nondim(m500_bothfree)
bothfree_array[4,1], bothfree_array[4,2], bothfree_array[4,3], bothfree_array[4,4], bothfree_array[4,5] = get_avg_forces_nondim(m1000_bothfree)

fixedOP_array = np.zeros((5,6))
fixedOP_array[0:5,0] = [50,250,375,500,675]
fixedOP_array[0,1], fixedOP_array[0,2], fixedOP_array[0,3], fixedOP_array[0,4], fixedOP_array[0,5] = get_avg_forces_nondim(m50_fixedOP)
fixedOP_array[1,1], fixedOP_array[1,2], fixedOP_array[1,3], fixedOP_array[1,4], fixedOP_array[1,5] = get_avg_forces_nondim(m250_fixedOP)
fixedOP_array[2,1], fixedOP_array[2,2], fixedOP_array[2,3], fixedOP_array[2,4], fixedOP_array[2,5] = get_avg_forces_nondim(m375_fixedOP)
fixedOP_array[3,1], fixedOP_array[3,2], fixedOP_array[3,3], fixedOP_array[3,4], fixedOP_array[3,5] = get_avg_forces_nondim(m500_fixedOP)
fixedOP_array[4,1], fixedOP_array[4,2], fixedOP_array[4,3], fixedOP_array[4,4], fixedOP_array[4,5] = get_avg_forces_nondim(m1000_fixedOP)


fig=plt.figure()
gs=GridSpec(3,2) 

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
for i in range(len(fixedSP_array)):

    # start with in-slab stresses
    plt.bar(fixedSP_array[i,0], (fixedSP_array[i,4]+fixedSP_array[i,5])/fixedSP_array[i,1], width=40, color='blue')
    
    if (fixedSP_array[i,4]+fixedSP_array[i,5])/fixedSP_array[i,1] < 0:
        bottom2 = 0
        thick2 = fixedSP_array[i,3]/fixedSP_array[i,1]
    else:
        bottom2 = (fixedSP_array[i,4]+fixedSP_array[i,5])/fixedSP_array[i,1]
        thick2 = ((fixedSP_array[i,4]+fixedSP_array[i,5])/fixedSP_array[i,1]) + fixedSP_array[i,3]/fixedSP_array[i,1] 

    plt.bar(fixedSP_array[i,0], thick2, bottom=bottom2, width=40, color='red')

plt.scatter(fixedSP_array[:,0],fixedSP_array[:,2]/fixedSP_array[:,1],s=20,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')

ax.set_xticks( [50,250, 375, 500,675] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='-',linewidth=0.75, zorder=0)
plt.axhline(y=1, color='gray',linestyle='--',linewidth=0.75, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("force/buoyancy",size=6)
# plt.xlim(-30,  1100); 
plt.ylim(-0.5,  1.5); 

ax=fig.add_subplot(gs[1,0])

# fixed SP models
for i in range(len(bothfree_array)):

    # start with in-slab stresses
    plt.bar(bothfree_array[i,0], (bothfree_array[i,4]+bothfree_array[i,5])/bothfree_array[i,1], width=40, color='blue')

    if (bothfree_array[i,4]+bothfree_array[i,5])/bothfree_array[i,1] < 0:
        bottom2 = 0
        thick2 = bothfree_array[i,3]/bothfree_array[i,1]
    else:
        bottom2 = (bothfree_array[i,4]+bothfree_array[i,5])/bothfree_array[i,1]
        thick2 =    bothfree_array[i,3]/bothfree_array[i,1]

    plt.bar(bothfree_array[i,0], thick2, bottom=bottom2, width=40, color='red')

plt.scatter(bothfree_array[:,0],bothfree_array[:,2]/bothfree_array[:,1],s=20,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')


ax.set_xticks( [50,250, 375, 500,675] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='-',linewidth=0.75, zorder=0)
plt.axhline(y=1, color='gray',linestyle='--',linewidth=0.75, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("force/buoyancy",size=6)
# plt.xlim(-30,  1100); 
plt.ylim(-0.5,  1.5); 


ax=fig.add_subplot(gs[2,0])

# fixed OP models
for i in range(len(fixedOP_array)):

    # start with in-slab stresses
    plt.bar(fixedOP_array[i,0], (fixedOP_array[i,4]+fixedOP_array[i,5])/fixedOP_array[i,1], width=40, color='blue')

    if (fixedOP_array[i,4]+fixedOP_array[i,5])/fixedOP_array[i,1] < 0:
        bottom2 = 0
        thick2 = fixedOP_array[i,3]/fixedOP_array[i,1]
    else:
        bottom2 = (fixedOP_array[i,4]+fixedOP_array[i,5])/fixedOP_array[i,1]
        thick2 = fixedOP_array[i,3]/fixedOP_array[i,1]

    if thick2 < 0:
        plt.bar(fixedOP_array[i,0], -1.0*thick2, bottom=(bottom2+thick2), width=33, color='red')
    else:
        plt.bar(fixedOP_array[i,0], thick2, bottom=bottom2, width=40, color='red')

plt.scatter(fixedOP_array[:,0],fixedOP_array[:,2]/fixedOP_array[:,1],s=20,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')


ax.set_xticks( [50,250, 375, 500,675] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='-',linewidth=0.75, zorder=0)
plt.axhline(y=1, color='gray',linestyle='--',linewidth=0.75, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("force/buoyancy",size=6)
# plt.xlim(-30,  1100); 
plt.ylim(-0.5,  1.5); 



# # 50
# dplot=30
# plt.scatter(50-dplot,mean_woshear_m50_fixedSP,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
# plt.plot([50-dplot,50-dplot], [mean_woshear_m50_fixedSP-stdev_woshear_m50_fixedSP,mean_woshear_m50_fixedSP+stdev_woshear_m50_fixedSP],  color='tan',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(50,mean_woshear_m50_bothfree,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([50,50], [mean_woshear_m50_bothfree-stdev_woshear_m50_bothfree,mean_woshear_m50_bothfree+stdev_woshear_m50_bothfree],  color='tan',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(50+dplot,mean_woshear_m50_fixedOP,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# plt.plot([50+dplot,50+dplot], [mean_woshear_m50_fixedOP-stdev_woshear_m50_fixedOP,mean_woshear_m50_fixedOP+stdev_woshear_m50_fixedOP],  color='tan',linewidth=3.5,zorder=2,alpha=0.5) 
# # 250
# plt.scatter(250-dplot,mean_woshear_m250_fixedSP,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([250-dplot,250-dplot], [mean_woshear_m250_fixedSP-stdev_woshear_m250_fixedSP,mean_woshear_m250_fixedSP+stdev_woshear_m250_fixedSP],  color='peru',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(250,mean_woshear_m250_bothfree,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([250,250], [mean_woshear_m250_bothfree-stdev_woshear_m250_bothfree,mean_woshear_m250_bothfree+stdev_woshear_m250_bothfree],  color='peru',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(250+dplot,mean_woshear_m250_fixedOP,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([250+dplot,250+dplot], [mean_woshear_m250_fixedOP-stdev_woshear_m250_fixedOP,mean_woshear_m250_fixedOP+stdev_woshear_m250_fixedOP],  color='peru',linewidth=3.5,zorder=2,alpha=0.5)
# # 375
# plt.scatter(375-dplot,mean_woshear_m375_fixedSP,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([375-dplot,375-dplot], [mean_woshear_m375_fixedSP-stdev_woshear_m375_fixedSP,mean_woshear_m375_fixedSP+stdev_woshear_m375_fixedSP],  color='firebrick',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(375,mean_woshear_m375_bothfree,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([375,375], [mean_woshear_m375_bothfree-stdev_woshear_m375_bothfree,mean_woshear_m375_bothfree+stdev_woshear_m375_bothfree],  color='firebrick',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(375+dplot,mean_woshear_m375_fixedOP,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([375+dplot,375+dplot], [mean_woshear_m375_fixedOP-stdev_woshear_m375_fixedOP,mean_woshear_m375_fixedOP+stdev_woshear_m375_fixedOP],  color='firebrick',linewidth=3.5,zorder=2,alpha=0.5)
# # 500
# plt.scatter(500-dplot,mean_woshear_m500_fixedSP,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([500-dplot,500-dplot], [mean_woshear_m500_fixedSP-stdev_woshear_m500_fixedSP,mean_woshear_m500_fixedSP+stdev_woshear_m500_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(500,mean_woshear_m500_bothfree,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([500,500], [mean_woshear_m500_bothfree-stdev_woshear_m500_bothfree,mean_woshear_m500_bothfree+stdev_woshear_m500_bothfree],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(500+dplot,mean_woshear_m500_fixedOP,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([500+dplot,500+dplot], [mean_woshear_m500_fixedOP-stdev_woshear_m500_fixedOP,mean_woshear_m500_fixedOP+stdev_woshear_m500_fixedOP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# # 1000
# plt.scatter(1000-dplot,mean_woshear_m1000_fixedSP,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([1000-dplot,1000-dplot], [mean_woshear_m1000_fixedSP-stdev_woshear_m1000_fixedSP,mean_woshear_m1000_fixedSP+stdev_woshear_m1000_fixedSP],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(1000,mean_woshear_m1000_bothfree,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([1000,1000], [mean_woshear_m1000_bothfree-stdev_woshear_m1000_bothfree,mean_woshear_m1000_bothfree+stdev_woshear_m1000_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(1000+dplot,mean_woshear_m1000_fixedOP,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([1000+dplot,1000+dplot], [mean_woshear_m1000_fixedOP-stdev_woshear_m1000_fixedOP,mean_woshear_m1000_fixedOP+stdev_woshear_m1000_fixedOP],  color='black',linewidth=3.5,zorder=2,alpha=0.5)


# ax.set_xticks( [50,250, 375, 500,1000] )
# ax.set_xticklabels( ['50','250','375','500','1000'] )
# plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
# plt.axhline(y=100, color='lightgray',linestyle='-',linewidth=0.75, zorder=0)
# ax.tick_params(axis='x', labelsize=6)
# ax.tick_params(axis='y', labelsize=6)
# plt.ylabel("misft (w/o shear)  [%]",size=6)
# plt.xlim(-30,  1100); 
# #plt.ylim(-8,  120); 
# plt.ylim(-20,  120); 

# fixed_aspect_ratio(0.8)

# #------

# ax=fig.add_subplot(gs[1,0])
# # 50
# plt.scatter(50-dplot,mean_wshear_m50_fixedSP,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
# plt.plot([50-dplot,50-dplot], [mean_wshear_m50_fixedSP-stdev_wshear_m50_fixedSP,mean_wshear_m50_fixedSP+stdev_wshear_m50_fixedSP],  color='tan',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(50,mean_wshear_m50_bothfree,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([50,50], [mean_wshear_m50_bothfree-stdev_wshear_m50_bothfree,mean_wshear_m50_bothfree+stdev_wshear_m50_bothfree],  color='tan',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(50+dplot,mean_wshear_m50_fixedOP,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# plt.plot([50+dplot,50+dplot], [mean_wshear_m50_fixedOP-stdev_wshear_m50_fixedOP,mean_wshear_m50_fixedOP+stdev_wshear_m50_fixedOP],  color='tan',linewidth=3.5,zorder=2,alpha=0.5) 
# # 250
# plt.scatter(250-dplot,mean_wshear_m250_fixedSP,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([250-dplot,250-dplot], [mean_wshear_m250_fixedSP-stdev_wshear_m250_fixedSP,mean_wshear_m250_fixedSP+stdev_wshear_m250_fixedSP],  color='peru',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(250,mean_wshear_m250_bothfree,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([250,250], [mean_wshear_m250_bothfree-stdev_wshear_m250_bothfree,mean_wshear_m250_bothfree+stdev_wshear_m250_bothfree],  color='peru',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(250+dplot,mean_wshear_m250_fixedOP,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([250+dplot,250+dplot], [mean_wshear_m250_fixedOP-stdev_wshear_m250_fixedOP,mean_wshear_m250_fixedOP+stdev_wshear_m250_fixedOP],  color='peru',linewidth=3.5,zorder=2,alpha=0.5)
# # 375
# plt.scatter(375-dplot,mean_wshear_m375_fixedSP,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([375-dplot,375-dplot], [mean_wshear_m375_fixedSP-stdev_wshear_m375_fixedSP,mean_wshear_m375_fixedSP+stdev_wshear_m375_fixedSP],  color='firebrick',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(375,mean_wshear_m375_bothfree,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([375,375], [mean_wshear_m375_bothfree-stdev_wshear_m375_bothfree,mean_wshear_m375_bothfree+stdev_wshear_m375_bothfree],  color='firebrick',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(375+dplot,mean_wshear_m375_fixedOP,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([375+dplot,375+dplot], [mean_wshear_m375_fixedOP-stdev_wshear_m375_fixedOP,mean_wshear_m375_fixedOP+stdev_wshear_m375_fixedOP],  color='firebrick',linewidth=3.5,zorder=2,alpha=0.5)
# # 500
# plt.scatter(500-dplot,mean_wshear_m500_fixedSP,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([500-dplot,500-dplot], [mean_wshear_m500_fixedSP-stdev_wshear_m500_fixedSP,mean_wshear_m500_fixedSP+stdev_wshear_m500_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(500,mean_wshear_m500_bothfree,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([500,500], [mean_wshear_m500_bothfree-stdev_wshear_m500_bothfree,mean_wshear_m500_bothfree+stdev_wshear_m500_bothfree],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(500+dplot,mean_wshear_m500_fixedOP,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([500+dplot,500+dplot], [mean_wshear_m500_fixedOP-stdev_wshear_m500_fixedOP,mean_wshear_m500_fixedOP+stdev_wshear_m500_fixedOP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# # 1000
# plt.scatter(1000-dplot,mean_wshear_m1000_fixedSP,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v') 
# plt.plot([1000-dplot,1000-dplot], [mean_wshear_m1000_fixedSP-stdev_wshear_m1000_fixedSP,mean_wshear_m1000_fixedSP+stdev_wshear_m1000_fixedSP],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(1000,mean_wshear_m1000_bothfree,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
# plt.plot([1000,1000], [mean_wshear_m1000_bothfree-stdev_wshear_m1000_bothfree,mean_wshear_m1000_bothfree+stdev_wshear_m1000_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(1000+dplot,mean_wshear_m1000_fixedOP,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^') 
# plt.plot([1000+dplot,1000+dplot], [mean_wshear_m1000_fixedOP-stdev_wshear_m1000_fixedOP,mean_wshear_m1000_fixedOP+stdev_wshear_m1000_fixedOP],  color='black',linewidth=3.5,zorder=2,alpha=0.5)


# ax.set_xticks( [50,250, 375, 500,1000] )
# ax.set_xticklabels( ['50','250','375','500','1000'] )
# plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
# ax.tick_params(axis='x', labelsize=6)
# ax.tick_params(axis='y', labelsize=6)
# plt.ylabel("misft (w/ shear)  [%]",size=6)
# plt.xlabel("slab strength",size=6)
# plt.xlim(-30,  1100); 
# plt.ylim(-20,  120); 
# plt.axhline(y=100, color='lightgray',linestyle='-',linewidth=0.75, zorder=0)
# fixed_aspect_ratio(0.8)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


