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
from functions import get_misfit_mean_and_stdev_thresh_time

analysis_depth  = float(sys.argv[1]) 
analysis_depth_dz = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

analysis_depth1 = 230e3
analysis_depth = 330e3
analysis_depth3 = 430e3

plot_name_png = ''.join(['plots/DP-comparisons/compilations/averaged-misfit_compilation_mature-subd.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/averaged-misfit_compilation_mature-subd.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])

# 50
name1_bothfree 	= "2D_compositional_subd_lower-res_new_50plates"
name1_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name1_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_50plates"
# 100
name2_bothfree 	="2D_compositional_subd_lower-res_new_WeakPlates"
name2_fixedSP    ="2D_compositional_subd_lower-res_new_FixedSP_WeakPlates"
name2_fixedOP    ="2D_compositional_subd_lower-res_new_FixedOP_WeakPlates"
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
# 2000
name6_bothfree ="2D_compositional_subd_lower-res_new_StiffPlates"
name6_fixedSP  ="2D_compositional_subd_lower-res_new_FixedSP_StiffPlates"
name6_fixedOP  ="2D_compositional_subd_lower-res_new_FixedOP_StiffPlates"

# text files
text50_bothfree 		= ''.join(['text_files/',name1_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedSP  		= ''.join(['text_files/',name1_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedOP  		= ''.join(['text_files/',name1_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text100_bothfree		= ''.join(['text_files/',name2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text100_fixedSP  		= ''.join(['text_files/',name2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text100_fixedOP  		= ''.join(['text_files/',name2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_bothfree    	= ''.join(['text_files/',name3_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedSP     	= ''.join(['text_files/',name3_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedOP		    = ''.join(['text_files/',name3_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_bothfree 		= ''.join(['text_files/',name4_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedSP  		= ''.join(['text_files/',name4_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedOP  		= ''.join(['text_files/',name4_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_bothfree		= ''.join(['text_files/',name5_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedSP  		= ''.join(['text_files/',name5_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedOP  		= ''.join(['text_files/',name5_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2000_bothfree    	= ''.join(['text_files/',name6_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2000_fixedSP     	= ''.join(['text_files/',name6_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2000_fixedOP	    = ''.join(['text_files/',name6_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


# load in models
m50_bothfree 	= np.loadtxt((text50_bothfree)) 
m50_fixedSP  	= np.loadtxt((text50_fixedSP))
m50_fixedOP  	= np.loadtxt((text50_fixedOP))
m100_bothfree 	= np.loadtxt((text100_bothfree)) 
m100_fixedSP  	= np.loadtxt((text100_fixedSP))
m100_fixedOP  	= np.loadtxt((text100_fixedOP))
m250_bothfree 	= np.loadtxt((text250_bothfree)) 
m250_fixedSP 	= np.loadtxt((text250_fixedSP))
m250_fixedOP 	= np.loadtxt((text250_fixedOP))
m500_bothfree 	= np.loadtxt((text500_bothfree)) 
m500_fixedSP  	= np.loadtxt((text500_fixedSP))
m500_fixedOP  	= np.loadtxt((text500_fixedOP))
m1000_bothfree 	= np.loadtxt((text1000_bothfree)) 
m1000_fixedSP  	= np.loadtxt((text1000_fixedSP))
m1000_fixedOP  	= np.loadtxt((text1000_fixedOP))
m2000_bothfree 	= np.loadtxt((text2000_bothfree)) 
m2000_fixedSP 	= np.loadtxt((text2000_fixedSP))
m2000_fixedOP 	= np.loadtxt((text2000_fixedOP))
# --

# get misfit means and st. deviations
thresold_t = 25.
mean_woshear_m50_bothfree, mean_wshear_m50_bothfree, stdev_woshear_m50_bothfree, stdev_wshear_m50_bothfree 	= get_misfit_mean_and_stdev_thresh_time(m50_bothfree,thresold_t)
mean_woshear_m50_fixedSP,  mean_wshear_m50_fixedSP,  stdev_woshear_m50_fixedSP,  stdev_wshear_m50_fixedSP 	= get_misfit_mean_and_stdev_thresh_time(m50_fixedSP,thresold_t)
mean_woshear_m50_fixedOP,  mean_wshear_m50_fixedOP,  stdev_woshear_m50_fixedOP,  stdev_wshear_m50_fixedOP 	= get_misfit_mean_and_stdev_thresh_time(m50_fixedOP,thresold_t)
mean_woshear_m100_bothfree, mean_wshear_m100_bothfree, stdev_woshear_m100_bothfree, stdev_wshear_m100_bothfree  = get_misfit_mean_and_stdev_thresh_time(m100_bothfree,thresold_t)
mean_woshear_m100_fixedSP,  mean_wshear_m100_fixedSP,  stdev_woshear_m100_fixedSP,  stdev_wshear_m100_fixedSP 	= get_misfit_mean_and_stdev_thresh_time(m100_fixedSP,thresold_t)
mean_woshear_m100_fixedOP,  mean_wshear_m100_fixedOP,  stdev_woshear_m100_fixedOP,  stdev_wshear_m100_fixedOP 	= get_misfit_mean_and_stdev_thresh_time(m100_fixedOP,thresold_t)
mean_woshear_m250_bothfree, mean_wshear_m250_bothfree, stdev_woshear_m250_bothfree, stdev_wshear_m250_bothfree 	= get_misfit_mean_and_stdev_thresh_time(m250_bothfree,thresold_t)
mean_woshear_m250_fixedSP,  mean_wshear_m250_fixedSP,  stdev_woshear_m250_fixedSP,  stdev_wshear_m250_fixedSP 	= get_misfit_mean_and_stdev_thresh_time(m250_fixedSP,thresold_t)
mean_woshear_m250_fixedOP,  mean_wshear_m250_fixedOP,  stdev_woshear_m250_fixedOP,  stdev_wshear_m250_fixedOP 	= get_misfit_mean_and_stdev_thresh_time(m250_fixedOP,thresold_t)
mean_woshear_m500_bothfree, mean_wshear_m500_bothfree, stdev_woshear_m500_bothfree, stdev_wshear_m500_bothfree 	= get_misfit_mean_and_stdev_thresh_time(m500_bothfree,thresold_t)
mean_woshear_m500_fixedSP,  mean_wshear_m500_fixedSP,  stdev_woshear_m500_fixedSP,  stdev_wshear_m500_fixedSP 	= get_misfit_mean_and_stdev_thresh_time(m500_fixedSP,thresold_t)
mean_woshear_m500_fixedOP,  mean_wshear_m500_fixedOP,  stdev_woshear_m500_fixedOP,  stdev_wshear_m500_fixedOP 	= get_misfit_mean_and_stdev_thresh_time(m500_fixedOP,thresold_t)
mean_woshear_m1000_bothfree, mean_wshear_m1000_bothfree, stdev_woshear_m1000_bothfree, stdev_wshear_m1000_bothfree 	= get_misfit_mean_and_stdev_thresh_time(m1000_bothfree,thresold_t)
mean_woshear_m1000_fixedSP,  mean_wshear_m1000_fixedSP,  stdev_woshear_m1000_fixedSP,  stdev_wshear_m1000_fixedSP 	= get_misfit_mean_and_stdev_thresh_time(m1000_fixedSP,thresold_t)
mean_woshear_m1000_fixedOP,  mean_wshear_m1000_fixedOP,  stdev_woshear_m1000_fixedOP,  stdev_wshear_m1000_fixedOP 	= get_misfit_mean_and_stdev_thresh_time(m1000_fixedOP,thresold_t)
mean_woshear_m2000_bothfree, mean_wshear_m2000_bothfree, stdev_woshear_m2000_bothfree, stdev_wshear_m2000_bothfree 	= get_misfit_mean_and_stdev_thresh_time(m2000_bothfree,thresold_t)
mean_woshear_m2000_fixedSP,  mean_wshear_m2000_fixedSP,  stdev_woshear_m2000_fixedSP,  stdev_wshear_m2000_fixedSP 	= get_misfit_mean_and_stdev_thresh_time(m2000_fixedSP,thresold_t)
mean_woshear_m2000_fixedOP,  mean_wshear_m2000_fixedOP,  stdev_woshear_m2000_fixedOP,  stdev_wshear_m2000_fixedOP 	= get_misfit_mean_and_stdev_thresh_time(m2000_fixedOP,thresold_t)


fig=plt.figure()
gs=GridSpec(2,1) 

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
# 50
dplot=47
plt.scatter(50-dplot,mean_woshear_m50_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.plot([50-dplot,50-dplot], [mean_woshear_m50_fixedSP-stdev_woshear_m50_fixedSP,mean_woshear_m50_fixedSP+stdev_woshear_m50_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(50,mean_woshear_m50_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([50,50], [mean_woshear_m50_bothfree-stdev_woshear_m50_bothfree,mean_woshear_m50_bothfree+stdev_woshear_m50_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(50+dplot,mean_woshear_m50_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3)
plt.plot([50+dplot,50+dplot], [mean_woshear_m50_fixedOP-stdev_woshear_m50_fixedOP,mean_woshear_m50_fixedOP+stdev_woshear_m50_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5) 
# # 100
# plt.scatter(100-dplot,mean_woshear_m100_fixedSP,s=22.5,color='maroon',zorder=3) 
# plt.plot([100-dplot,100-dplot], [mean_woshear_m100_fixedSP-stdev_woshear_m100_fixedSP,mean_woshear_m100_fixedSP+stdev_woshear_m100_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(100,mean_woshear_m100_bothfree,s=22.5,color='black',zorder=3)
# plt.plot([100,100], [mean_woshear_m100_bothfree-stdev_woshear_m100_bothfree,mean_woshear_m100_bothfree+stdev_woshear_m100_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5) 
# plt.scatter(100+dplot,mean_woshear_m100_fixedOP,s=22.5,color='slategray',zorder=3) 
# plt.plot([100+dplot,100+dplot], [mean_woshear_m100_fixedOP-stdev_woshear_m100_fixedOP,mean_woshear_m100_fixedOP+stdev_woshear_m100_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# 250
plt.scatter(250-dplot,mean_woshear_m250_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([250-dplot,250-dplot], [mean_woshear_m250_fixedSP-stdev_woshear_m250_fixedSP,mean_woshear_m250_fixedSP+stdev_woshear_m250_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(250,mean_woshear_m250_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([250,250], [mean_woshear_m250_bothfree-stdev_woshear_m250_bothfree,mean_woshear_m250_bothfree+stdev_woshear_m250_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(250+dplot,mean_woshear_m250_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([250+dplot,250+dplot], [mean_woshear_m250_fixedOP-stdev_woshear_m250_fixedOP,mean_woshear_m250_fixedOP+stdev_woshear_m250_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# 500
plt.scatter(500-dplot,mean_woshear_m500_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([500-dplot,500-dplot], [mean_woshear_m500_fixedSP-stdev_woshear_m500_fixedSP,mean_woshear_m500_fixedSP+stdev_woshear_m500_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(500,mean_woshear_m500_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([500,500], [mean_woshear_m500_bothfree-stdev_woshear_m500_bothfree,mean_woshear_m500_bothfree+stdev_woshear_m500_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(500+dplot,mean_woshear_m500_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([500+dplot,500+dplot], [mean_woshear_m500_fixedOP-stdev_woshear_m500_fixedOP,mean_woshear_m500_fixedOP+stdev_woshear_m500_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# 1000
plt.scatter(1000-dplot,mean_woshear_m1000_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([1000-dplot,1000-dplot], [mean_woshear_m1000_fixedSP-stdev_woshear_m1000_fixedSP,mean_woshear_m1000_fixedSP+stdev_woshear_m1000_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(1000,mean_woshear_m1000_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([1000,1000], [mean_woshear_m1000_bothfree-stdev_woshear_m1000_bothfree,mean_woshear_m1000_bothfree+stdev_woshear_m1000_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(1000+dplot,mean_woshear_m1000_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([1000+dplot,1000+dplot], [mean_woshear_m1000_fixedOP-stdev_woshear_m1000_fixedOP,mean_woshear_m1000_fixedOP+stdev_woshear_m1000_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)

ax.set_xticks( [50,250,500,1000] )
ax.set_xticklabels( ['50','250','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("misft (w/o shear)  [MPa]",size=6)
plt.xlim(-30,  1100); 
plt.ylim(-1,  15); 
fixed_aspect_ratio(0.8)

#------

ax=fig.add_subplot(gs[1,0])
# 50
plt.scatter(50-dplot,mean_wshear_m50_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.plot([50-dplot,50-dplot], [mean_wshear_m50_fixedSP-stdev_wshear_m50_fixedSP,mean_wshear_m50_fixedSP+stdev_wshear_m50_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(50,mean_wshear_m50_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([50,50], [mean_wshear_m50_bothfree-stdev_wshear_m50_bothfree,mean_wshear_m50_bothfree+stdev_wshear_m50_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(50+dplot,mean_wshear_m50_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3)
plt.plot([50+dplot,50+dplot], [mean_wshear_m50_fixedOP-stdev_wshear_m50_fixedOP,mean_wshear_m50_fixedOP+stdev_wshear_m50_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5) 
# # 100
# plt.scatter(100-dplot,mean_wshear_m100_fixedSP,s=22.5,color='maroon',zorder=3) 
# plt.plot([100-dplot,100-dplot], [mean_wshear_m100_fixedSP-stdev_wshear_m100_fixedSP,mean_wshear_m100_fixedSP+stdev_wshear_m100_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
# plt.scatter(100,mean_wshear_m100_bothfree,s=22.5,color='black',zorder=3)
# plt.plot([100,100], [mean_wshear_m100_bothfree-stdev_wshear_m100_bothfree,mean_wshear_m100_bothfree+stdev_wshear_m100_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5) 
# plt.scatter(100+dplot,mean_wshear_m100_fixedOP,s=22.5,color='slategray',zorder=3) 
# plt.plot([100+dplot,100+dplot], [mean_wshear_m100_fixedOP-stdev_wshear_m100_fixedOP,mean_wshear_m100_fixedOP+stdev_wshear_m100_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# 250
plt.scatter(250-dplot,mean_wshear_m250_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([250-dplot,250-dplot], [mean_wshear_m250_fixedSP-stdev_wshear_m250_fixedSP,mean_wshear_m250_fixedSP+stdev_wshear_m250_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(250,mean_wshear_m250_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([250,250], [mean_wshear_m250_bothfree-stdev_wshear_m250_bothfree,mean_wshear_m250_bothfree+stdev_wshear_m250_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(250+dplot,mean_wshear_m250_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([250+dplot,250+dplot], [mean_wshear_m250_fixedOP-stdev_wshear_m250_fixedOP,mean_wshear_m250_fixedOP+stdev_wshear_m250_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# 500
plt.scatter(500-dplot,mean_wshear_m500_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([500-dplot,500-dplot], [mean_wshear_m500_fixedSP-stdev_wshear_m500_fixedSP,mean_wshear_m500_fixedSP+stdev_wshear_m500_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(500,mean_wshear_m500_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([500,500], [mean_wshear_m500_bothfree-stdev_wshear_m500_bothfree,mean_wshear_m500_bothfree+stdev_wshear_m500_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(500+dplot,mean_wshear_m500_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([500+dplot,500+dplot], [mean_wshear_m500_fixedOP-stdev_wshear_m500_fixedOP,mean_wshear_m500_fixedOP+stdev_wshear_m500_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# 1000
plt.scatter(1000-dplot,mean_wshear_m1000_fixedSP,s=22.5,color='maroon',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([1000-dplot,1000-dplot], [mean_wshear_m1000_fixedSP-stdev_wshear_m1000_fixedSP,mean_wshear_m1000_fixedSP+stdev_wshear_m1000_fixedSP],  color='maroon',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(1000,mean_wshear_m1000_bothfree,s=22.5,color='black',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([1000,1000], [mean_wshear_m1000_bothfree-stdev_wshear_m1000_bothfree,mean_wshear_m1000_bothfree+stdev_wshear_m1000_bothfree],  color='black',linewidth=3.5,zorder=2,alpha=0.5)
plt.scatter(1000+dplot,mean_wshear_m1000_fixedOP,s=22.5,color='slategray',edgecolor='black',linewidth=0.25,zorder=3) 
plt.plot([1000+dplot,1000+dplot], [mean_wshear_m1000_fixedOP-stdev_wshear_m1000_fixedOP,mean_wshear_m1000_fixedOP+stdev_wshear_m1000_fixedOP],  color='slategray',linewidth=3.5,zorder=2,alpha=0.5)
# # 2000

ax.set_xticks( [50,250,500,1000] )
ax.set_xticklabels( ['50','250','500','1000'] )
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("misft (w/ shear)  [MPa]",size=6)
plt.xlabel("slab strength",size=6)
plt.xlim(-30,  1100); 
plt.ylim(-1,  15); 
fixed_aspect_ratio(0.8)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


