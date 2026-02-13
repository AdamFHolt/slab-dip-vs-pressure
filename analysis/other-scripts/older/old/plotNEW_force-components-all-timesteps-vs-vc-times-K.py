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
m250_fixedSP 		= np.loadtxt((text250_fixedSP))
m250_fixedOP 		= np.loadtxt((text250_fixedOP))
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

K_ind    = 11
dK_ind   = 12
ss_ind   = 17
sn_ind   = 6
anal_ind = 4
DP_ind   = 3
dip_ind  = 5
vc_ind   = 20

plot_name_png = ''.join(['plots/DP-comparisons/compilations/NEW/force-components-all-vs-vc-times-K.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/NEW/force-components-all-vs-vc-times-K.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])
x_label = "vc * K [cm/km yr]"



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



# plot model DP vs. x-axis variable
fig=plt.figure()

ax=fig.add_subplot(gs[0,0])
# 50
plt.scatter(m50_fixedSP[:,K_ind]*m50_fixedSP[:,vc_ind],m50_fixedSP[:,DP_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[:,K_ind]*m50_bothfree[:,vc_ind],m50_bothfree[:,DP_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[:,K_ind]*m50_fixedOP[:,vc_ind],m50_fixedOP[:,DP_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[:,K_ind]*m250_fixedSP[:,vc_ind],m250_fixedSP[:,DP_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[:,K_ind]*m250_bothfree[:,vc_ind],m250_bothfree[:,DP_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[:,K_ind]*m250_fixedOP[:,vc_ind],m250_fixedOP[:,DP_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[:,K_ind]*m375_fixedSP[:,vc_ind],m375_fixedSP[:,DP_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[:,K_ind]*m375_bothfree[:,vc_ind],m375_bothfree[:,DP_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[:,K_ind]*m375_fixedOP[:,vc_ind],m375_fixedOP[:,DP_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[:,K_ind]*m500_fixedSP[:,vc_ind],m500_fixedSP[:,DP_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[:,K_ind]*m500_bothfree[:,vc_ind],m500_bothfree[:,DP_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[:,K_ind]*m500_fixedOP[:,vc_ind],m500_fixedOP[:,DP_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[:,K_ind]*m1000_fixedSP[:,vc_ind],m1000_fixedSP[:,DP_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[:,K_ind]*m1000_bothfree[:,vc_ind],m1000_bothfree[:,DP_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[:,K_ind]*m1000_fixedOP[:,vc_ind],m1000_fixedOP[:,DP_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("model DP  [MPa]",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-5,  30); 

# plot analytical DP vs. x-axis variable

ax=fig.add_subplot(gs[0,1])

# 50
plt.scatter(m50_fixedSP[:,K_ind]*m50_fixedSP[:,vc_ind],m50_fixedSP[:,anal_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[:,K_ind]*m50_bothfree[:,vc_ind],m50_bothfree[:,anal_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[:,K_ind]*m50_fixedOP[:,vc_ind],m50_fixedOP[:,anal_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[:,K_ind]*m250_fixedSP[:,vc_ind],m250_fixedSP[:,anal_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[:,K_ind]*m250_bothfree[:,vc_ind],m250_bothfree[:,anal_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[:,K_ind]*m250_fixedOP[:,vc_ind],m250_fixedOP[:,anal_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[:,K_ind]*m375_fixedSP[:,vc_ind],m375_fixedSP[:,anal_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[:,K_ind]*m375_bothfree[:,vc_ind],m375_bothfree[:,anal_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[:,K_ind]*m375_fixedOP[:,vc_ind],m375_fixedOP[:,anal_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[:,K_ind]*m500_fixedSP[:,vc_ind],m500_fixedSP[:,anal_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[:,K_ind]*m500_bothfree[:,vc_ind],m500_bothfree[:,anal_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[:,K_ind]*m500_fixedOP[:,vc_ind],m500_fixedOP[:,anal_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[:,K_ind]*m1000_fixedSP[:,vc_ind],m1000_fixedSP[:,anal_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[:,K_ind]*m1000_bothfree[:,vc_ind],m1000_bothfree[:,anal_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[:,K_ind]*m1000_fixedOP[:,vc_ind],m1000_fixedOP[:,anal_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("analytical DP  [MPa]",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-5,  30); 


# plot model slab shear stress vs. x-axis variable

ax=fig.add_subplot(gs[1,0])

# 50
plt.scatter(m50_fixedSP[:,K_ind]*m50_fixedSP[:,vc_ind],-1.0*m50_fixedSP[:,sn_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[:,K_ind]*m50_bothfree[:,vc_ind],-1.0*m50_bothfree[:,sn_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[:,K_ind]*m50_fixedOP[:,vc_ind],-1.0*m50_fixedOP[:,sn_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[:,K_ind]*m250_fixedSP[:,vc_ind],-1.0*m250_fixedSP[:,sn_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[:,K_ind]*m250_bothfree[:,vc_ind],-1.0*m250_bothfree[:,sn_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[:,K_ind]*m250_fixedOP[:,vc_ind],-1.0*m250_fixedOP[:,sn_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[:,K_ind]*m375_fixedSP[:,vc_ind],-1.0*m375_fixedSP[:,sn_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[:,K_ind]*m375_bothfree[:,vc_ind],-1.0*m375_bothfree[:,sn_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[:,K_ind]*m375_fixedOP[:,vc_ind],-1.0*m375_fixedOP[:,sn_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[:,K_ind]*m500_fixedSP[:,vc_ind],-1.0*m500_fixedSP[:,sn_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[:,K_ind]*m500_bothfree[:,vc_ind],-1.0*m500_bothfree[:,sn_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[:,K_ind]*m500_fixedOP[:,vc_ind],-1.0*m500_fixedOP[:,sn_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[:,K_ind]*m1000_fixedSP[:,vc_ind],-1.0*m1000_fixedSP[:,sn_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[:,K_ind]*m1000_bothfree[:,vc_ind],-1.0*m1000_bothfree[:,sn_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[:,K_ind]*m1000_fixedOP[:,vc_ind],-1.0*m1000_fixedOP[:,sn_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("-1 * s-shear [MPa]",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-15,  25); 



# plot model slab normal stress vs. x-axis variable

ax=fig.add_subplot(gs[1,1])

# 50
plt.scatter(m50_fixedSP[:,K_ind]*m50_fixedSP[:,vc_ind],-1.0*m50_fixedSP[:,ss_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[:,K_ind]*m50_bothfree[:,vc_ind],-1.0*m50_bothfree[:,ss_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[:,K_ind]*m50_fixedOP[:,vc_ind],-1.0*m50_fixedOP[:,ss_ind]/1.e6,s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[:,K_ind]*m250_fixedSP[:,vc_ind],-1.0*m250_fixedSP[:,ss_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[:,K_ind]*m250_bothfree[:,vc_ind],-1.0*m250_bothfree[:,ss_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[:,K_ind]*m250_fixedOP[:,vc_ind],-1.0*m250_fixedOP[:,ss_ind]/1.e6,s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[:,K_ind]*m375_fixedSP[:,vc_ind],-1.0*m375_fixedSP[:,ss_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[:,K_ind]*m375_bothfree[:,vc_ind],-1.0*m375_bothfree[:,ss_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[:,K_ind]*m375_fixedOP[:,vc_ind],-1.0*m375_fixedOP[:,ss_ind]/1.e6,s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[:,K_ind]*m500_fixedSP[:,vc_ind],-1.0*m500_fixedSP[:,ss_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[:,K_ind]*m500_bothfree[:,vc_ind],-1.0*m500_bothfree[:,ss_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[:,K_ind]*m500_fixedOP[:,vc_ind],-1.0*m500_fixedOP[:,ss_ind]/1.e6,s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[:,K_ind]*m1000_fixedSP[:,vc_ind],-1.0*m1000_fixedSP[:,ss_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[:,K_ind]*m1000_bothfree[:,vc_ind],-1.0*m1000_bothfree[:,ss_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[:,K_ind]*m1000_fixedOP[:,vc_ind],-1.0*m1000_fixedOP[:,ss_ind]/1.e6,s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("-1 * s-norm [MPa]",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-15,  25); 

# plot model DP misfit vs. x-axis variable

ax=fig.add_subplot(gs[2,0])

# 50
plt.scatter(m50_fixedSP[:,K_ind]*m50_fixedSP[:,vc_ind],100.*((m50_fixedSP[:,DP_ind]-m50_fixedSP[:,anal_ind])/(m50_fixedSP[:,anal_ind])),s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[:,K_ind]*m50_bothfree[:,vc_ind],100.*((m50_bothfree[:,DP_ind]-m50_bothfree[:,anal_ind])/(m50_bothfree[:,anal_ind])),s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[:,K_ind]*m50_fixedOP[:,vc_ind],100.*((m50_fixedOP[:,DP_ind]-m50_fixedOP[:,anal_ind])/(m50_fixedOP[:,anal_ind])),s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[:,K_ind]*m250_fixedSP[:,vc_ind],100.*((m250_fixedSP[:,DP_ind]-m250_fixedSP[:,anal_ind])/(m250_fixedSP[:,anal_ind])),s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[:,K_ind]*m250_bothfree[:,vc_ind],100.*((m250_bothfree[:,DP_ind]-m250_bothfree[:,anal_ind])/(m250_bothfree[:,anal_ind])),s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[:,K_ind]*m250_fixedOP[:,vc_ind],100.*((m250_fixedOP[:,DP_ind]-m250_fixedOP[:,anal_ind])/(m250_fixedOP[:,anal_ind])),s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[:,K_ind]*m375_fixedSP[:,vc_ind],100.*((m375_fixedSP[:,DP_ind]-m375_fixedSP[:,anal_ind])/(m375_fixedSP[:,anal_ind])),s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[:,K_ind]*m375_bothfree[:,vc_ind],100.*((m375_bothfree[:,DP_ind]-m375_bothfree[:,anal_ind])/(m375_bothfree[:,anal_ind])),s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[:,K_ind]*m375_fixedOP[:,vc_ind],100.*((m375_fixedOP[:,DP_ind]-m375_fixedOP[:,anal_ind])/(m375_fixedOP[:,anal_ind])),s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[:,K_ind]*m500_fixedSP[:,vc_ind],100.*((m500_fixedSP[:,DP_ind]-m500_fixedSP[:,anal_ind])/(m500_fixedSP[:,anal_ind])),s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[:,K_ind]*m500_bothfree[:,vc_ind],100.*((m500_bothfree[:,DP_ind]-m500_bothfree[:,anal_ind])/(m500_bothfree[:,anal_ind])),s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[:,K_ind]*m500_fixedOP[:,vc_ind],100.*((m500_fixedOP[:,DP_ind]-m500_fixedOP[:,anal_ind])/(m500_fixedOP[:,anal_ind])),s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[:,K_ind]*m1000_fixedSP[:,vc_ind],100.*((m1000_fixedSP[:,DP_ind]-m1000_fixedSP[:,anal_ind])/(m1000_fixedSP[:,anal_ind])),s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[:,K_ind]*m1000_bothfree[:,vc_ind],100.*((m1000_bothfree[:,DP_ind]-m1000_bothfree[:,anal_ind])/(m1000_bothfree[:,anal_ind])),s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[:,K_ind]*m1000_fixedOP[:,vc_ind],100.*((m1000_fixedOP[:,DP_ind]-m1000_fixedOP[:,anal_ind])/(m1000_fixedOP[:,anal_ind])),s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("DP misfit [%]",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-100,  100); 


# plot model full stress misfit vs. x-axis variable

ax=fig.add_subplot(gs[2,1])

# 50
plt.scatter(m50_fixedSP[:,K_ind]*m50_fixedSP[:,vc_ind],100.*((m50_fixedSP[:,DP_ind]-m50_fixedSP[:,anal_ind]-m50_fixedSP[:,ss_ind]-m50_fixedSP[:,sn_ind])/(m50_fixedSP[:,anal_ind])),s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[:,K_ind]*m50_bothfree[:,vc_ind],100.*((m50_bothfree[:,DP_ind]-m50_bothfree[:,anal_ind]-m50_bothfree[:,ss_ind]-m50_bothfree[:,sn_ind])/(m50_bothfree[:,anal_ind])),s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[:,K_ind]*m50_fixedOP[:,vc_ind],100.*((m50_fixedOP[:,DP_ind]-m50_fixedOP[:,anal_ind]-m50_fixedOP[:,ss_ind]-m50_fixedOP[:,sn_ind])/(m50_fixedOP[:,anal_ind])),s=20,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[:,K_ind]*m250_fixedSP[:,vc_ind],100.*((m250_fixedSP[:,DP_ind]-m250_fixedSP[:,anal_ind]-m250_fixedSP[:,ss_ind]-m250_fixedSP[:,sn_ind])/(m250_fixedSP[:,anal_ind])),s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[:,K_ind]*m250_bothfree[:,vc_ind],100.*((m250_bothfree[:,DP_ind]-m250_bothfree[:,anal_ind]-m250_bothfree[:,ss_ind]-m250_bothfree[:,sn_ind])/(m250_bothfree[:,anal_ind])),s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[:,K_ind]*m250_fixedOP[:,vc_ind],100.*((m250_fixedOP[:,DP_ind]-m250_fixedOP[:,anal_ind]-m250_fixedOP[:,ss_ind]-m250_fixedOP[:,sn_ind])/(m250_fixedOP[:,anal_ind])),s=20,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[:,K_ind]*m375_fixedSP[:,vc_ind],100.*((m375_fixedSP[:,DP_ind]-m375_fixedSP[:,anal_ind]-m375_fixedSP[:,ss_ind]-m375_fixedSP[:,sn_ind])/(m375_fixedSP[:,anal_ind])),s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[:,K_ind]*m375_bothfree[:,vc_ind],100.*((m375_bothfree[:,DP_ind]-m375_bothfree[:,anal_ind]-m375_bothfree[:,ss_ind]-m375_bothfree[:,sn_ind])/(m375_bothfree[:,anal_ind])),s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[:,K_ind]*m375_fixedOP[:,vc_ind],100.*((m375_fixedOP[:,DP_ind]-m375_fixedOP[:,anal_ind]-m375_fixedOP[:,ss_ind]-m375_fixedOP[:,sn_ind])/(m375_fixedOP[:,anal_ind])),s=20,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[:,K_ind]*m500_fixedSP[:,vc_ind],100.*((m500_fixedSP[:,DP_ind]-m500_fixedSP[:,anal_ind]-m500_fixedSP[:,ss_ind]-m500_fixedSP[:,sn_ind])/(m500_fixedSP[:,anal_ind])),s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[:,K_ind]*m500_bothfree[:,vc_ind],100.*((m500_bothfree[:,DP_ind]-m500_bothfree[:,anal_ind]-m500_bothfree[:,ss_ind]-m500_bothfree[:,sn_ind])/(m500_bothfree[:,anal_ind])),s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[:,K_ind]*m500_fixedOP[:,vc_ind],100.*((m500_fixedOP[:,DP_ind]-m500_fixedOP[:,anal_ind]-m500_fixedOP[:,ss_ind]-m500_fixedOP[:,sn_ind])/(m500_fixedOP[:,anal_ind])),s=20,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[:,K_ind]*m1000_fixedSP[:,vc_ind],100.*((m1000_fixedSP[:,DP_ind]-m1000_fixedSP[:,anal_ind]-m1000_fixedSP[:,ss_ind]-m1000_fixedSP[:,sn_ind])/(m1000_fixedSP[:,anal_ind])),s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[:,K_ind]*m1000_bothfree[:,vc_ind],100.*((m1000_bothfree[:,DP_ind]-m1000_bothfree[:,anal_ind]-m1000_bothfree[:,ss_ind]-m1000_bothfree[:,sn_ind])/(m1000_bothfree[:,anal_ind])),s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[:,K_ind]*m1000_fixedOP[:,vc_ind],100.*((m1000_fixedOP[:,DP_ind]-m1000_fixedOP[:,anal_ind]-m1000_fixedOP[:,ss_ind]-m1000_fixedOP[:,sn_ind])/(m1000_fixedOP[:,anal_ind])),s=20,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("full stress misfit [%]",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-100,  100); 

plt.subplots_adjust(wspace=0.5)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')
