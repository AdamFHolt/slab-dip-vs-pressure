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
from functions import get_misfit_mean_and_stdev, get_curvature_mean_and_stdev, get_misfit_mean_and_stdev_nondim

analysis_depth  = float(sys.argv[1]) 
analysis_depth_dz = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

tactual_min = 14 # first time step to use
tmin = tactual_min - 8

plot_name_png = ''.join(['plots/DP-comparisons/compilations/averaged/misfits-vs-viscosity.K-colored.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/averaged/misfits-vs-viscosity.K-colored.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.pdf'])

name1_bothfree 	= "2D_compositional_subd_lower-res_new_50plates"
name1_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name1_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_50plates"
name3_bothfree 	= "2D_compositional_subd_lower-res_new_250plates"
name3_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_250plates"
name3_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_250plates"
name4_bothfree 	= "2D_compositional_subd_lower-res_new2"
name4_fixedSP  	= "2D_compositional_subd_FixedSP_lower-res_new2"
name4_fixedOP  	= "2D_compositional_subd_FixedOP_lower-res_new"
name5_bothfree = "2D_compositional_subd_lower-res_new_1000plates"
name5_fixedSP  = "2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
name5_fixedOP  = "2D_compositional_subd_lower-res_new_FixedOP_1000plates"
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


# get misfit means and st. deviations
mean_woshear_m50_bothfree, mean_wshear_m50_bothfree, stdev_woshear_m50_bothfree, stdev_wshear_m50_bothfree 	= get_misfit_mean_and_stdev(m50_bothfree,tmin)
mean_woshear_m50_fixedSP,  mean_wshear_m50_fixedSP,  stdev_woshear_m50_fixedSP,  stdev_wshear_m50_fixedSP 	= get_misfit_mean_and_stdev(m50_fixedSP,tmin)
mean_woshear_m50_fixedOP,  mean_wshear_m50_fixedOP,  stdev_woshear_m50_fixedOP,  stdev_wshear_m50_fixedOP 	= get_misfit_mean_and_stdev(m50_fixedOP,tmin)
mean_woshear_m250_bothfree, mean_wshear_m250_bothfree, stdev_woshear_m250_bothfree, stdev_wshear_m250_bothfree 	= get_misfit_mean_and_stdev(m250_bothfree,tmin)
mean_woshear_m250_fixedSP,  mean_wshear_m250_fixedSP,  stdev_woshear_m250_fixedSP,  stdev_wshear_m250_fixedSP 	= get_misfit_mean_and_stdev(m250_fixedSP,tmin)
mean_woshear_m250_fixedOP,  mean_wshear_m250_fixedOP,  stdev_woshear_m250_fixedOP,  stdev_wshear_m250_fixedOP 	= get_misfit_mean_and_stdev(m250_fixedOP,tmin)
mean_woshear_m500_bothfree, mean_wshear_m500_bothfree, stdev_woshear_m500_bothfree, stdev_wshear_m500_bothfree 	= get_misfit_mean_and_stdev(m500_bothfree,tmin)
mean_woshear_m500_fixedSP,  mean_wshear_m500_fixedSP,  stdev_woshear_m500_fixedSP,  stdev_wshear_m500_fixedSP 	= get_misfit_mean_and_stdev(m500_fixedSP,tmin)
mean_woshear_m500_fixedOP,  mean_wshear_m500_fixedOP,  stdev_woshear_m500_fixedOP,  stdev_wshear_m500_fixedOP 	= get_misfit_mean_and_stdev(m500_fixedOP,tmin)
mean_woshear_m1000_bothfree, mean_wshear_m1000_bothfree, stdev_woshear_m1000_bothfree, stdev_wshear_m1000_bothfree 	= get_misfit_mean_and_stdev(m1000_bothfree,tmin)
mean_woshear_m1000_fixedSP,  mean_wshear_m1000_fixedSP,  stdev_woshear_m1000_fixedSP,  stdev_wshear_m1000_fixedSP 	= get_misfit_mean_and_stdev(m1000_fixedSP,tmin)
mean_woshear_m1000_fixedOP,  mean_wshear_m1000_fixedOP,  stdev_woshear_m1000_fixedOP,  stdev_wshear_m1000_fixedOP 	= get_misfit_mean_and_stdev(m1000_fixedOP,tmin)
mean_woshear_m375_bothfree, mean_wshear_m375_bothfree, stdev_woshear_m375_bothfree, stdev_wshear_m375_bothfree 	= get_misfit_mean_and_stdev(m375_bothfree,tmin)
mean_woshear_m375_fixedSP,  mean_wshear_m375_fixedSP,  stdev_woshear_m375_fixedSP,  stdev_wshear_m375_fixedSP 	= get_misfit_mean_and_stdev(m375_fixedSP,tmin)
mean_woshear_m375_fixedOP,  mean_wshear_m375_fixedOP,  stdev_woshear_m375_fixedOP,  stdev_wshear_m375_fixedOP 	= get_misfit_mean_and_stdev(m375_fixedOP,tmin)

# get curvature means and st. deviations
mean_m50_bothfree,   stdev_m50_bothfree   = get_curvature_mean_and_stdev(m50_bothfree,tmin)
mean_m50_fixedSP,    stdev_m50_fixedSP    = get_curvature_mean_and_stdev(m50_fixedSP,tmin)
mean_m50_fixedOP,    stdev_m50_fixedOP    = get_curvature_mean_and_stdev(m50_fixedOP,tmin)
mean_m250_bothfree,  stdev_m250_bothfree  = get_curvature_mean_and_stdev(m250_bothfree,tmin)
mean_m250_fixedSP,   stdev_m250_fixedSP   = get_curvature_mean_and_stdev(m250_fixedSP,tmin)
mean_m250_fixedOP,   stdev_m250_fixedOP   = get_curvature_mean_and_stdev(m250_fixedOP,tmin)
mean_m500_bothfree,  stdev_m500_bothfree  = get_curvature_mean_and_stdev(m500_bothfree,tmin)
mean_m500_fixedSP,   stdev_m500_fixedSP   = get_curvature_mean_and_stdev(m500_fixedSP,tmin)
mean_m500_fixedOP,   stdev_m500_fixedOP   = get_curvature_mean_and_stdev(m500_fixedOP,tmin)
mean_m1000_bothfree, stdev_m1000_bothfree = get_curvature_mean_and_stdev(m1000_bothfree,tmin)
mean_m1000_fixedSP,  stdev_m1000_fixedSP  = get_curvature_mean_and_stdev(m1000_fixedSP,tmin)
mean_m1000_fixedOP,  stdev_m1000_fixedOP  = get_curvature_mean_and_stdev(m1000_fixedOP,tmin)
mean_m375_bothfree, stdev_m375_bothfree = get_curvature_mean_and_stdev(m375_bothfree,tmin)
mean_m375_fixedSP,  stdev_m375_fixedSP  = get_curvature_mean_and_stdev(m375_fixedSP,tmin)
mean_m375_fixedOP,  stdev_m375_fixedOP  = get_curvature_mean_and_stdev(m375_fixedOP,tmin)

fig=plt.figure()
gs=GridSpec(2,2) 
dplot=30



color_map = cm.get_cmap('Greys')
color_min = 0
#color_max = 0.0008
color_max = 0.001
interval  = color_max/20.
norm = matplotlib.colors.BoundaryNorm(np.arange(color_min,color_max+interval,interval), color_map.N)

ax=fig.add_subplot(gs[0,0])

#---

im = plt.scatter(50-dplot,mean_woshear_m50_fixedSP,s=20,c=mean_m50_fixedSP,cmap=color_map,norm=norm,edgecolor='tan',linewidth=0.5,zorder=3,marker='v')
plt.scatter(50,mean_woshear_m50_bothfree,s=20,c=mean_m50_bothfree,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3)
plt.scatter(50+dplot,mean_woshear_m50_fixedOP,s=20,c=mean_m50_fixedOP,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3,marker='^')

plt.plot([50-dplot,50-dplot], [mean_woshear_m50_fixedSP-stdev_woshear_m50_fixedSP,mean_woshear_m50_fixedSP+stdev_woshear_m50_fixedSP],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50,50], [mean_woshear_m50_bothfree-stdev_woshear_m50_bothfree,mean_woshear_m50_bothfree+stdev_woshear_m50_bothfree],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50+dplot,50+dplot], [mean_woshear_m50_fixedOP-stdev_woshear_m50_fixedOP,mean_woshear_m50_fixedOP+stdev_woshear_m50_fixedOP],  color='tan',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(250-dplot,mean_woshear_m250_fixedSP,s=20,c=mean_m250_fixedSP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='v')
plt.scatter(250,mean_woshear_m250_bothfree,s=20,c=mean_m250_bothfree,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3)
plt.scatter(250+dplot,mean_woshear_m250_fixedOP,s=20,c=mean_m250_fixedOP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='^')

plt.plot([250-dplot,250-dplot], [mean_woshear_m250_fixedSP-stdev_woshear_m250_fixedSP,mean_woshear_m250_fixedSP+stdev_woshear_m250_fixedSP],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250,250], [mean_woshear_m250_bothfree-stdev_woshear_m250_bothfree,mean_woshear_m250_bothfree+stdev_woshear_m250_bothfree],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250+dplot,250+dplot], [mean_woshear_m250_fixedOP-stdev_woshear_m250_fixedOP,mean_woshear_m250_fixedOP+stdev_woshear_m250_fixedOP],  color='peru',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(375-dplot,mean_woshear_m375_fixedSP,s=20,c=mean_m375_fixedSP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='v')
plt.scatter(375,mean_woshear_m375_bothfree,s=20,c=mean_m375_bothfree,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3)
plt.scatter(375+dplot,mean_woshear_m375_fixedOP,s=20,c=mean_m375_fixedOP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='^')

plt.plot([375-dplot,375-dplot], [mean_woshear_m375_fixedSP-stdev_woshear_m375_fixedSP,mean_woshear_m375_fixedSP+stdev_woshear_m375_fixedSP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375,375], [mean_woshear_m375_bothfree-stdev_woshear_m375_bothfree,mean_woshear_m375_bothfree+stdev_woshear_m375_bothfree],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375+dplot,375+dplot], [mean_woshear_m375_fixedOP-stdev_woshear_m375_fixedOP,mean_woshear_m375_fixedOP+stdev_woshear_m375_fixedOP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(500-dplot,mean_woshear_m500_fixedSP,s=20,c=mean_m500_fixedSP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='v')
plt.scatter(500,mean_woshear_m500_bothfree,s=20,c=mean_m500_bothfree,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3)
plt.scatter(500+dplot,mean_woshear_m500_fixedOP,s=20,c=mean_m500_fixedOP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='^')

plt.plot([500-dplot,500-dplot], [mean_woshear_m500_fixedSP-stdev_woshear_m500_fixedSP,mean_woshear_m500_fixedSP+stdev_woshear_m500_fixedSP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500,500], [mean_woshear_m500_bothfree-stdev_woshear_m500_bothfree,mean_woshear_m500_bothfree+stdev_woshear_m500_bothfree],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500+dplot,500+dplot], [mean_woshear_m500_fixedOP-stdev_woshear_m500_fixedOP,mean_woshear_m500_fixedOP+stdev_woshear_m500_fixedOP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(1000-dplot,mean_woshear_m1000_fixedSP,s=20,c=mean_m1000_fixedSP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='v')
plt.scatter(1000,mean_woshear_m1000_bothfree,s=20,c=mean_m1000_bothfree,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3)
plt.scatter(1000+dplot,mean_woshear_m1000_fixedOP,s=20,c=mean_m1000_fixedOP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='^')

plt.plot([1000-dplot,1000-dplot], [mean_woshear_m1000_fixedSP-stdev_woshear_m1000_fixedSP,mean_woshear_m1000_fixedSP+stdev_woshear_m1000_fixedSP],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000,1000], [mean_woshear_m1000_bothfree-stdev_woshear_m1000_bothfree,mean_woshear_m1000_bothfree+stdev_woshear_m1000_bothfree],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000+dplot,1000+dplot], [mean_woshear_m1000_fixedOP-stdev_woshear_m1000_fixedOP,mean_woshear_m1000_fixedOP+stdev_woshear_m1000_fixedOP],  color='black',linewidth=2,zorder=2,alpha=0.5)

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("misfit (w/o slab stress)  [MPa]",size=6)
plt.ylim(-2,  15); 
plt.xlim(-30,  1100);
fixed_aspect_ratio(1)
ax.set_xticks( [50,250, 375, 500,1000] )
ax.set_xticklabels( ['50','250','375','500','1000'] )

cbar = plt.colorbar(im, cax = fig.add_axes([0.19, 0.82, 0.08, 0.008]), orientation='horizontal', ticks=[0, color_max/2, color_max],extend='max') # left, base, length left-right, length top-bott
cbar.cmap.set_over('k')
cbar.ax.tick_params(axis='x',labelsize=4.5,pad=1)
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.xaxis.set_label_position('bottom')
cbar.ax.set_title('K  [1/km]',size=6)

# with shear 

ax=fig.add_subplot(gs[0,1])


plt.scatter(50-dplot,mean_wshear_m50_fixedSP,s=20,c=mean_m50_fixedSP,cmap=color_map,norm=norm,edgecolor='tan',linewidth=0.5,zorder=3,marker='v')
plt.scatter(50,mean_wshear_m50_bothfree,s=20,c=mean_m50_bothfree,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3)
plt.scatter(50+dplot,mean_wshear_m50_fixedOP,s=20,c=mean_m50_fixedOP,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3,marker='^')

plt.plot([50-dplot,50-dplot], [mean_wshear_m50_fixedSP-stdev_wshear_m50_fixedSP,mean_wshear_m50_fixedSP+stdev_wshear_m50_fixedSP],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50,50], [mean_wshear_m50_bothfree-stdev_wshear_m50_bothfree,mean_wshear_m50_bothfree+stdev_wshear_m50_bothfree],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50+dplot,50+dplot], [mean_wshear_m50_fixedOP-stdev_wshear_m50_fixedOP,mean_wshear_m50_fixedOP+stdev_wshear_m50_fixedOP],  color='tan',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(250-dplot,mean_wshear_m250_fixedSP,s=20,c=mean_m250_fixedSP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='v')
plt.scatter(250,mean_wshear_m250_bothfree,s=20,c=mean_m250_bothfree,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3)
plt.scatter(250+dplot,mean_wshear_m250_fixedOP,s=20,c=mean_m250_fixedOP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='^')

plt.plot([250-dplot,250-dplot], [mean_wshear_m250_fixedSP-stdev_wshear_m250_fixedSP,mean_wshear_m250_fixedSP+stdev_wshear_m250_fixedSP],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250,250], [mean_wshear_m250_bothfree-stdev_wshear_m250_bothfree,mean_wshear_m250_bothfree+stdev_wshear_m250_bothfree],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250+dplot,250+dplot], [mean_wshear_m250_fixedOP-stdev_wshear_m250_fixedOP,mean_wshear_m250_fixedOP+stdev_wshear_m250_fixedOP],  color='peru',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(375-dplot,mean_wshear_m375_fixedSP,s=20,c=mean_m375_fixedSP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='v')
plt.scatter(375,mean_wshear_m375_bothfree,s=20,c=mean_m375_bothfree,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3)
plt.scatter(375+dplot,mean_wshear_m375_fixedOP,s=20,c=mean_m375_fixedOP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='^')

plt.plot([375-dplot,375-dplot], [mean_wshear_m375_fixedSP-stdev_wshear_m375_fixedSP,mean_wshear_m375_fixedSP+stdev_wshear_m375_fixedSP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375,375], [mean_wshear_m375_bothfree-stdev_wshear_m375_bothfree,mean_wshear_m375_bothfree+stdev_wshear_m375_bothfree],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375+dplot,375+dplot], [mean_wshear_m375_fixedOP-stdev_wshear_m375_fixedOP,mean_wshear_m375_fixedOP+stdev_wshear_m375_fixedOP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(500-dplot,mean_wshear_m500_fixedSP,s=20,c=mean_m500_fixedSP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='v')
plt.scatter(500,mean_wshear_m500_bothfree,s=20,c=mean_m500_bothfree,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3)
plt.scatter(500+dplot,mean_wshear_m500_fixedOP,s=20,c=mean_m500_fixedOP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='^')

plt.plot([500-dplot,500-dplot], [mean_wshear_m500_fixedSP-stdev_wshear_m500_fixedSP,mean_wshear_m500_fixedSP+stdev_wshear_m500_fixedSP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500,500], [mean_wshear_m500_bothfree-stdev_wshear_m500_bothfree,mean_wshear_m500_bothfree+stdev_wshear_m500_bothfree],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500+dplot,500+dplot], [mean_wshear_m500_fixedOP-stdev_wshear_m500_fixedOP,mean_wshear_m500_fixedOP+stdev_wshear_m500_fixedOP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(1000-dplot,mean_wshear_m1000_fixedSP,s=20,c=mean_m1000_fixedSP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='v')
plt.scatter(1000,mean_wshear_m1000_bothfree,s=20,c=mean_m1000_bothfree,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3)
plt.scatter(1000+dplot,mean_wshear_m1000_fixedOP,s=20,c=mean_m1000_fixedOP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='^')

plt.plot([1000-dplot,1000-dplot], [mean_wshear_m1000_fixedSP-stdev_wshear_m1000_fixedSP,mean_wshear_m1000_fixedSP+stdev_wshear_m1000_fixedSP],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000,1000], [mean_wshear_m1000_bothfree-stdev_wshear_m1000_bothfree,mean_wshear_m1000_bothfree+stdev_wshear_m1000_bothfree],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000+dplot,1000+dplot], [mean_wshear_m1000_fixedOP-stdev_wshear_m1000_fixedOP,mean_wshear_m1000_fixedOP+stdev_wshear_m1000_fixedOP],  color='black',linewidth=2,zorder=2,alpha=0.5)


    
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("misfit (w/  slab stress)  [MPa]",size=6)
plt.ylim(-2,  15); 
plt.xlim(-30,  1100);
fixed_aspect_ratio(1)
ax.set_xticks( [50,250, 375, 500,1000] )
ax.set_xticklabels( ['50','250','375','500','1000'] )

# now non-dim
mean_woshear_m50_bothfree, mean_wshear_m50_bothfree, stdev_woshear_m50_bothfree, stdev_wshear_m50_bothfree  = get_misfit_mean_and_stdev_nondim(m50_bothfree,tmin)
mean_woshear_m50_fixedSP,  mean_wshear_m50_fixedSP,  stdev_woshear_m50_fixedSP,  stdev_wshear_m50_fixedSP   = get_misfit_mean_and_stdev_nondim(m50_fixedSP,tmin)
mean_woshear_m50_fixedOP,  mean_wshear_m50_fixedOP,  stdev_woshear_m50_fixedOP,  stdev_wshear_m50_fixedOP   = get_misfit_mean_and_stdev_nondim(m50_fixedOP,tmin)
mean_woshear_m250_both
free, mean_wshear_m250_bothfree, stdev_woshear_m250_bothfree, stdev_wshear_m250_bothfree  = get_misfit_mean_and_stdev_nondim(m250_bothfree,tmin)
mean_woshear_m250_fixedSP,  mean_wshear_m250_fixedSP,  stdev_woshear_m250_fixedSP,  stdev_wshear_m250_fixedSP   = get_misfit_mean_and_stdev_nondim(m250_fixedSP,tmin)
mean_woshear_m250_fixedOP,  mean_wshear_m250_fixedOP,  stdev_woshear_m250_fixedOP,  stdev_wshear_m250_fixedOP   = get_misfit_mean_and_stdev_nondim(m250_fixedOP,tmin)
mean_woshear_m500_bothfree, mean_wshear_m500_bothfree, stdev_woshear_m500_bothfree, stdev_wshear_m500_bothfree  = get_misfit_mean_and_stdev_nondim(m500_bothfree,tmin)
mean_woshear_m500_fixedSP,  mean_wshear_m500_fixedSP,  stdev_woshear_m500_fixedSP,  stdev_wshear_m500_fixedSP   = get_misfit_mean_and_stdev_nondim(m500_fixedSP,tmin)
mean_woshear_m500_fixedOP,  mean_wshear_m500_fixedOP,  stdev_woshear_m500_fixedOP,  stdev_wshear_m500_fixedOP   = get_misfit_mean_and_stdev_nondim(m500_fixedOP,tmin)
mean_woshear_m1000_bothfree, mean_wshear_m1000_bothfree, stdev_woshear_m1000_bothfree, stdev_wshear_m1000_bothfree  = get_misfit_mean_and_stdev_nondim(m1000_bothfree,tmin)
mean_woshear_m1000_fixedSP,  mean_wshear_m1000_fixedSP,  stdev_woshear_m1000_fixedSP,  stdev_wshear_m1000_fixedSP   = get_misfit_mean_and_stdev_nondim(m1000_fixedSP,tmin)
mean_woshear_m1000_fixedOP,  mean_wshear_m1000_fixedOP,  stdev_woshear_m1000_fixedOP,  stdev_wshear_m1000_fixedOP   = get_misfit_mean_and_stdev_nondim(m1000_fixedOP,tmin)
mean_woshear_m375_bothfree, mean_wshear_m375_bothfree, stdev_woshear_m375_bothfree, stdev_wshear_m375_bothfree  = get_misfit_mean_and_stdev_nondim(m375_bothfree,tmin)
mean_woshear_m375_fixedSP,  mean_wshear_m375_fixedSP,  stdev_woshear_m375_fixedSP,  stdev_wshear_m375_fixedSP   = get_misfit_mean_and_stdev_nondim(m375_fixedSP,tmin)
mean_woshear_m375_fixedOP,  mean_wshear_m375_fixedOP,  stdev_woshear_m375_fixedOP,  stdev_wshear_m375_fixedOP   = get_misfit_mean_and_stdev_nondim(m375_fixedOP,tmin)


ax=fig.add_subplot(gs[1,0])

#---

im = plt.scatter(50-dplot,mean_woshear_m50_fixedSP,s=20,c=mean_m50_fixedSP,cmap=color_map,norm=norm,edgecolor='tan',linewidth=0.5,zorder=3,marker='v')
plt.scatter(50,mean_woshear_m50_bothfree,s=20,c=mean_m50_bothfree,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3)
plt.scatter(50+dplot,mean_woshear_m50_fixedOP,s=20,c=mean_m50_fixedOP,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3,marker='^')

plt.plot([50-dplot,50-dplot], [mean_woshear_m50_fixedSP-stdev_woshear_m50_fixedSP,mean_woshear_m50_fixedSP+stdev_woshear_m50_fixedSP],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50,50], [mean_woshear_m50_bothfree-stdev_woshear_m50_bothfree,mean_woshear_m50_bothfree+stdev_woshear_m50_bothfree],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50+dplot,50+dplot], [mean_woshear_m50_fixedOP-stdev_woshear_m50_fixedOP,mean_woshear_m50_fixedOP+stdev_woshear_m50_fixedOP],  color='tan',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(250-dplot,mean_woshear_m250_fixedSP,s=20,c=mean_m250_fixedSP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='v')
plt.scatter(250,mean_woshear_m250_bothfree,s=20,c=mean_m250_bothfree,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3)
plt.scatter(250+dplot,mean_woshear_m250_fixedOP,s=20,c=mean_m250_fixedOP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='^')

plt.plot([250-dplot,250-dplot], [mean_woshear_m250_fixedSP-stdev_woshear_m250_fixedSP,mean_woshear_m250_fixedSP+stdev_woshear_m250_fixedSP],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250,250], [mean_woshear_m250_bothfree-stdev_woshear_m250_bothfree,mean_woshear_m250_bothfree+stdev_woshear_m250_bothfree],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250+dplot,250+dplot], [mean_woshear_m250_fixedOP-stdev_woshear_m250_fixedOP,mean_woshear_m250_fixedOP+stdev_woshear_m250_fixedOP],  color='peru',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(375-dplot,mean_woshear_m375_fixedSP,s=20,c=mean_m375_fixedSP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='v')
plt.scatter(375,mean_woshear_m375_bothfree,s=20,c=mean_m375_bothfree,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3)
plt.scatter(375+dplot,mean_woshear_m375_fixedOP,s=20,c=mean_m375_fixedOP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='^')

plt.plot([375-dplot,375-dplot], [mean_woshear_m375_fixedSP-stdev_woshear_m375_fixedSP,mean_woshear_m375_fixedSP+stdev_woshear_m375_fixedSP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375,375], [mean_woshear_m375_bothfree-stdev_woshear_m375_bothfree,mean_woshear_m375_bothfree+stdev_woshear_m375_bothfree],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375+dplot,375+dplot], [mean_woshear_m375_fixedOP-stdev_woshear_m375_fixedOP,mean_woshear_m375_fixedOP+stdev_woshear_m375_fixedOP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(500-dplot,mean_woshear_m500_fixedSP,s=20,c=mean_m500_fixedSP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='v')
plt.scatter(500,mean_woshear_m500_bothfree,s=20,c=mean_m500_bothfree,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3)
plt.scatter(500+dplot,mean_woshear_m500_fixedOP,s=20,c=mean_m500_fixedOP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='^')

plt.plot([500-dplot,500-dplot], [mean_woshear_m500_fixedSP-stdev_woshear_m500_fixedSP,mean_woshear_m500_fixedSP+stdev_woshear_m500_fixedSP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500,500], [mean_woshear_m500_bothfree-stdev_woshear_m500_bothfree,mean_woshear_m500_bothfree+stdev_woshear_m500_bothfree],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500+dplot,500+dplot], [mean_woshear_m500_fixedOP-stdev_woshear_m500_fixedOP,mean_woshear_m500_fixedOP+stdev_woshear_m500_fixedOP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(1000-dplot,mean_woshear_m1000_fixedSP,s=20,c=mean_m1000_fixedSP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='v')
plt.scatter(1000,mean_woshear_m1000_bothfree,s=20,c=mean_m1000_bothfree,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3)
plt.scatter(1000+dplot,mean_woshear_m1000_fixedOP,s=20,c=mean_m1000_fixedOP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='^')

plt.plot([1000-dplot,1000-dplot], [mean_woshear_m1000_fixedSP-stdev_woshear_m1000_fixedSP,mean_woshear_m1000_fixedSP+stdev_woshear_m1000_fixedSP],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000,1000], [mean_woshear_m1000_bothfree-stdev_woshear_m1000_bothfree,mean_woshear_m1000_bothfree+stdev_woshear_m1000_bothfree],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000+dplot,1000+dplot], [mean_woshear_m1000_fixedOP-stdev_woshear_m1000_fixedOP,mean_woshear_m1000_fixedOP+stdev_woshear_m1000_fixedOP],  color='black',linewidth=2,zorder=2,alpha=0.5)

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("misfit (w/o slab stress)  [%]",size=6)
plt.xlabel("viscosity  [Pa s]",size=6)
plt.ylim(-8,  120);
plt.xlim(-30,  1100);
fixed_aspect_ratio(1)
ax.set_xticks( [50,250, 375, 500,1000] )
ax.set_xticklabels( ['50','250','375','500','1000'] )


# with shear 
ax=fig.add_subplot(gs[1,1])

plt.scatter(50-dplot,mean_wshear_m50_fixedSP,s=20,c=mean_m50_fixedSP,cmap=color_map,norm=norm,edgecolor='tan',linewidth=0.5,zorder=3,marker='v')
plt.scatter(50,mean_wshear_m50_bothfree,s=20,c=mean_m50_bothfree,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3)
plt.scatter(50+dplot,mean_wshear_m50_fixedOP,s=20,c=mean_m50_fixedOP,cmap=color_map,norm=norm, edgecolor='tan',linewidth=0.5,zorder=3,marker='^')

plt.plot([50-dplot,50-dplot], [mean_wshear_m50_fixedSP-stdev_wshear_m50_fixedSP,mean_wshear_m50_fixedSP+stdev_wshear_m50_fixedSP],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50,50], [mean_wshear_m50_bothfree-stdev_wshear_m50_bothfree,mean_wshear_m50_bothfree+stdev_wshear_m50_bothfree],  color='tan',linewidth=2,zorder=2,alpha=0.5)
plt.plot([50+dplot,50+dplot], [mean_wshear_m50_fixedOP-stdev_wshear_m50_fixedOP,mean_wshear_m50_fixedOP+stdev_wshear_m50_fixedOP],  color='tan',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(250-dplot,mean_wshear_m250_fixedSP,s=20,c=mean_m250_fixedSP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='v')
plt.scatter(250,mean_wshear_m250_bothfree,s=20,c=mean_m250_bothfree,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3)
plt.scatter(250+dplot,mean_wshear_m250_fixedOP,s=20,c=mean_m250_fixedOP,cmap=color_map,norm=norm,edgecolor='peru',linewidth=0.5,zorder=3,marker='^')

plt.plot([250-dplot,250-dplot], [mean_wshear_m250_fixedSP-stdev_wshear_m250_fixedSP,mean_wshear_m250_fixedSP+stdev_wshear_m250_fixedSP],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250,250], [mean_wshear_m250_bothfree-stdev_wshear_m250_bothfree,mean_wshear_m250_bothfree+stdev_wshear_m250_bothfree],  color='peru',linewidth=2,zorder=2,alpha=0.5)
plt.plot([250+dplot,250+dplot], [mean_wshear_m250_fixedOP-stdev_wshear_m250_fixedOP,mean_wshear_m250_fixedOP+stdev_wshear_m250_fixedOP],  color='peru',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(375-dplot,mean_wshear_m375_fixedSP,s=20,c=mean_m375_fixedSP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='v')
plt.scatter(375,mean_wshear_m375_bothfree,s=20,c=mean_m375_bothfree,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3)
plt.scatter(375+dplot,mean_wshear_m375_fixedOP,s=20,c=mean_m375_fixedOP,cmap=color_map,norm=norm,edgecolor='firebrick',linewidth=0.5,zorder=3,marker='^')

plt.plot([375-dplot,375-dplot], [mean_wshear_m375_fixedSP-stdev_wshear_m375_fixedSP,mean_wshear_m375_fixedSP+stdev_wshear_m375_fixedSP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375,375], [mean_wshear_m375_bothfree-stdev_wshear_m375_bothfree,mean_wshear_m375_bothfree+stdev_wshear_m375_bothfree],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)
plt.plot([375+dplot,375+dplot], [mean_wshear_m375_fixedOP-stdev_wshear_m375_fixedOP,mean_wshear_m375_fixedOP+stdev_wshear_m375_fixedOP],  color='firebrick',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(500-dplot,mean_wshear_m500_fixedSP,s=20,c=mean_m500_fixedSP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='v')
plt.scatter(500,mean_wshear_m500_bothfree,s=20,c=mean_m500_bothfree,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3)
plt.scatter(500+dplot,mean_wshear_m500_fixedOP,s=20,c=mean_m500_fixedOP,cmap=color_map,norm=norm,edgecolor='maroon',linewidth=0.5,zorder=3,marker='^')

plt.plot([500-dplot,500-dplot], [mean_wshear_m500_fixedSP-stdev_wshear_m500_fixedSP,mean_wshear_m500_fixedSP+stdev_wshear_m500_fixedSP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500,500], [mean_wshear_m500_bothfree-stdev_wshear_m500_bothfree,mean_wshear_m500_bothfree+stdev_wshear_m500_bothfree],  color='maroon',linewidth=2,zorder=2,alpha=0.5)
plt.plot([500+dplot,500+dplot], [mean_wshear_m500_fixedOP-stdev_wshear_m500_fixedOP,mean_wshear_m500_fixedOP+stdev_wshear_m500_fixedOP],  color='maroon',linewidth=2,zorder=2,alpha=0.5)

#---

plt.scatter(1000-dplot,mean_wshear_m1000_fixedSP,s=20,c=mean_m1000_fixedSP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='v')
plt.scatter(1000,mean_wshear_m1000_bothfree,s=20,c=mean_m1000_bothfree,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3)
plt.scatter(1000+dplot,mean_wshear_m1000_fixedOP,s=20,c=mean_m1000_fixedOP,cmap=color_map,norm=norm,edgecolor='black',linewidth=0.5,zorder=3,marker='^')

plt.plot([1000-dplot,1000-dplot], [mean_wshear_m1000_fixedSP-stdev_wshear_m1000_fixedSP,mean_wshear_m1000_fixedSP+stdev_wshear_m1000_fixedSP],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000,1000], [mean_wshear_m1000_bothfree-stdev_wshear_m1000_bothfree,mean_wshear_m1000_bothfree+stdev_wshear_m1000_bothfree],  color='black',linewidth=2,zorder=2,alpha=0.5)
plt.plot([1000+dplot,1000+dplot], [mean_wshear_m1000_fixedOP-stdev_wshear_m1000_fixedOP,mean_wshear_m1000_fixedOP+stdev_wshear_m1000_fixedOP],  color='black',linewidth=2,zorder=2,alpha=0.5)


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("misfit (w/  slab stress)  [%]",size=6)
plt.xlabel("viscosity  [Pa s]",size=6)
plt.ylim(-8,  120);
plt.xlim(-30,  1100);
fixed_aspect_ratio(1)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
#plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')
