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

tactual_min = 11 # first time step to use
tmin = tactual_min - 8

plot_name_png = ''.join(['plots/DP-comparisons/compilations/force-components-all-vs-scaling.non-dim.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/force-components-all-vs-scaling.non-dim.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.pdf'])

# 50
name1_bothfree 	= "2D_compositional_subd_lower-res_new_50plates"
name1_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name1_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_50plates"
# 250
name3_bothfree 	= "2D_compositional_subd_lower-res_new_250plates"
name3_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_250plates"
name3_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_250plates"
# 500
name4_bothfree 	= "2D_compositional_subd_lower-res_new2"
name4_fixedSP  	= "2D_compositional_subd_FixedSP_lower-res_new2"
name4_fixedOP  	= "2D_compositional_subd_FixedOP_lower-res_new"
# 1000
name5_bothfree = "2D_compositional_subd_lower-res_new_1000plates"
name5_fixedSP  = "2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
name5_fixedOP  = "2D_compositional_subd_lower-res_new_FixedOP_1000plates"
# 375
name7_bothfree 	= "2D_compositional_subd_lower-res_new_375plates"
name7_fixedSP  	= "2D_compositional_subd_lower-res_new_FixedSP_375plates"
name7_fixedOP  	= "2D_compositional_subd_lower-res_new_FixedOP_375plates"


# text files
text50_bothfree 		= ''.join(['text_files/',name1_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedSP  		= ''.join(['text_files/',name1_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedOP  		= ''.join(['text_files/',name1_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_bothfree    	= ''.join(['text_files/',name3_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedSP     	= ''.join(['text_files/',name3_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedOP		    = ''.join(['text_files/',name3_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_bothfree 		= ''.join(['text_files/',name4_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedSP  		= ''.join(['text_files/',name4_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedOP  		= ''.join(['text_files/',name4_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_bothfree		= ''.join(['text_files/',name5_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedSP  		= ''.join(['text_files/',name5_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedOP  		= ''.join(['text_files/',name5_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_bothfree    	= ''.join(['text_files/',name7_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedSP     	= ''.join(['text_files/',name7_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedOP	    	= ''.join(['text_files/',name7_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

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

K_ind    = 11
ss_ind   = 17
sn_ind   = 6
anal_ind = 4
DP_ind   = 3
dip_ind  = 5
vc_ind   = 19

x_ind = K_ind
x_label = "V K visc [Pa]"
m50_bothfree[:,K_ind] 	=  m50_bothfree[:,K_ind]   * 1.e3 * 50 * m50_bothfree[:,vc_ind]
m50_fixedSP[:,K_ind] 	=  m50_fixedSP[:,K_ind]    * 1.e3 * 50 * m50_fixedSP[:,vc_ind]
m50_fixedOP[:,K_ind] 	=  m50_fixedOP[:,K_ind]    * 1.e3 * 50 * m50_fixedOP[:,vc_ind]
m250_bothfree[:,K_ind] 	=  m250_bothfree[:,K_ind]  * 1.e3 * 250 * m250_bothfree[:,vc_ind]
m250_fixedSP[:,K_ind] 	=  m250_fixedSP[:,K_ind]   * 1.e3 * 250 * m250_fixedSP[:,vc_ind]
m250_fixedOP[:,K_ind] 	=  m250_fixedOP[:,K_ind]   * 1.e3 * 250 * m250_fixedOP[:,vc_ind]
m375_bothfree[:,K_ind] 	=  m375_bothfree[:,K_ind]  * 1.e3 * 375 * m375_bothfree[:,vc_ind]
m375_fixedSP[:,K_ind] 	=  m375_fixedSP[:,K_ind]   * 1.e3 * 375 * m375_fixedSP[:,vc_ind]
m375_fixedOP[:,K_ind] 	=  m375_fixedOP[:,K_ind]   * 1.e3 * 375 * m375_fixedOP[:,vc_ind]
m500_bothfree[:,K_ind] 	=  m500_bothfree[:,K_ind]  * 1.e3 * 500 * m500_bothfree[:,vc_ind]
m500_fixedSP[:,K_ind] 	=  m500_fixedSP[:,K_ind]   * 1.e3 * 500 * m500_fixedSP[:,vc_ind]
m500_fixedOP[:,K_ind] 	=  m500_fixedOP[:,K_ind]   * 1.e3 * 500 * m500_fixedOP[:,vc_ind]
m1000_bothfree[:,K_ind] =  m1000_bothfree[:,K_ind] * 1.e3 * 1000 * m1000_bothfree[:,vc_ind]
m1000_fixedSP[:,K_ind] 	=  m1000_fixedSP[:,K_ind]  * 1.e3 * 1000 * m1000_fixedSP[:,vc_ind]
m1000_fixedOP[:,K_ind] 	=  m1000_fixedOP[:,K_ind]  * 1.e3 * 1000 * m1000_fixedOP[:,vc_ind]

# m50_bothfree[:,K_ind] 	=  m50_bothfree[:,K_ind]   * 1.e3 * 50 * m50_bothfree[:,vc_ind]      * (1/np.sin(np.deg2rad(m50_bothfree[:,dip_ind])))
# m50_fixedSP[:,K_ind] 	=  m50_fixedSP[:,K_ind]    * 1.e3 * 50 * m50_fixedSP[:,vc_ind]       * (1/np.sin(np.deg2rad(m50_fixedSP[:,dip_ind])))
# m50_fixedOP[:,K_ind] 	=  m50_fixedOP[:,K_ind]    * 1.e3 * 50 * m50_fixedOP[:,vc_ind]       * (1/np.sin(np.deg2rad(m50_fixedOP[:,dip_ind])))
# m250_bothfree[:,K_ind] 	=  m250_bothfree[:,K_ind]  * 1.e3 * 250 * m250_bothfree[:,vc_ind]    * (1/np.sin(np.deg2rad(m250_bothfree[:,dip_ind])))
# m250_fixedSP[:,K_ind] 	=  m250_fixedSP[:,K_ind]   * 1.e3 * 250 * m250_fixedSP[:,vc_ind]     * (1/np.sin(np.deg2rad(m250_fixedSP[:,dip_ind])))
# m250_fixedOP[:,K_ind] 	=  m250_fixedOP[:,K_ind]   * 1.e3 * 250 * m250_fixedOP[:,vc_ind]     * (1/np.sin(np.deg2rad(m250_fixedOP[:,dip_ind])))
# m375_bothfree[:,K_ind] 	=  m375_bothfree[:,K_ind]  * 1.e3 * 375 * m375_bothfree[:,vc_ind]    * (1/np.sin(np.deg2rad(m375_bothfree[:,dip_ind])))
# m375_fixedSP[:,K_ind] 	=  m375_fixedSP[:,K_ind]   * 1.e3 * 375 * m375_fixedSP[:,vc_ind]     * (1/np.sin(np.deg2rad(m375_fixedSP[:,dip_ind])))
# m375_fixedOP[:,K_ind] 	=  m375_fixedOP[:,K_ind]   * 1.e3 * 375 * m375_fixedOP[:,vc_ind]     * (1/np.sin(np.deg2rad(m375_fixedOP[:,dip_ind])))
# m500_bothfree[:,K_ind] 	=  m500_bothfree[:,K_ind]  * 1.e3 * 500 * m500_bothfree[:,vc_ind]    * (1/np.sin(np.deg2rad(m500_bothfree[:,dip_ind])))
# m500_fixedSP[:,K_ind] 	=  m500_fixedSP[:,K_ind]   * 1.e3 * 500 * m500_fixedSP[:,vc_ind]     * (1/np.sin(np.deg2rad(m500_fixedSP[:,dip_ind])))
# m500_fixedOP[:,K_ind] 	=  m500_fixedOP[:,K_ind]   * 1.e3 * 500 * m500_fixedOP[:,vc_ind]     * (1/np.sin(np.deg2rad(m500_fixedOP[:,dip_ind])))
# m1000_bothfree[:,K_ind] =  m1000_bothfree[:,K_ind] * 1.e3 * 1000 * m1000_bothfree[:,vc_ind]  * (1/np.sin(np.deg2rad(m1000_bothfree[:,dip_ind])))
# m1000_fixedSP[:,K_ind] 	=  m1000_fixedSP[:,K_ind]  * 1.e3 * 1000 * m1000_fixedSP[:,vc_ind]   * (1/np.sin(np.deg2rad(m1000_fixedSP[:,dip_ind])))
# m1000_fixedOP[:,K_ind] 	=  m1000_fixedOP[:,K_ind]  * 1.e3 * 1000 * m1000_fixedOP[:,vc_ind]   * (1/np.sin(np.deg2rad(m1000_fixedOP[:,dip_ind])))


gs=GridSpec(2,2) 

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

fig=plt.figure()

# plot model full stress misfit vs. x-axis variable

ax=fig.add_subplot(gs[0,0])

# 50
plt.scatter(m50_fixedSP[tmin:,x_ind], (m50_fixedSP[tmin:,DP_ind]-m50_fixedSP[tmin:,anal_ind]+m50_fixedSP[tmin:,ss_ind]-m50_fixedSP[tmin:,sn_ind])/m50_fixedSP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[tmin:,x_ind], (m50_bothfree[tmin:,DP_ind]-m50_bothfree[tmin:,anal_ind]+m50_bothfree[tmin:,ss_ind]-m50_bothfree[tmin:,sn_ind])/m50_bothfree[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[tmin:,x_ind], (m50_fixedOP[tmin:,DP_ind]-m50_fixedOP[tmin:,anal_ind]+m50_fixedOP[tmin:,ss_ind]-m50_fixedOP[tmin:,sn_ind])/m50_fixedOP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[tmin:,x_ind], (m250_fixedSP[tmin:,DP_ind]-m250_fixedSP[tmin:,anal_ind]+m250_fixedSP[tmin:,ss_ind]-m250_fixedSP[tmin:,sn_ind])/m250_fixedSP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[tmin:,x_ind], (m250_bothfree[tmin:,DP_ind]-m250_bothfree[tmin:,anal_ind]+m250_bothfree[tmin:,ss_ind]-m250_bothfree[tmin:,sn_ind])/m250_bothfree[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[tmin:,x_ind], (m250_fixedOP[tmin:,DP_ind]-m250_fixedOP[tmin:,anal_ind]+m250_fixedOP[tmin:,ss_ind]-m250_fixedOP[tmin:,sn_ind])/m250_fixedOP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[tmin:,x_ind], (m375_fixedSP[tmin:,DP_ind]-m375_fixedSP[tmin:,anal_ind]+m375_fixedSP[tmin:,ss_ind]-m375_fixedSP[tmin:,sn_ind])/m375_fixedSP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[tmin:,x_ind], (m375_bothfree[tmin:,DP_ind]-m375_bothfree[tmin:,anal_ind]+m375_bothfree[tmin:,ss_ind]-m375_bothfree[tmin:,sn_ind])/m375_bothfree[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[tmin:,x_ind], (m375_fixedOP[tmin:,DP_ind]-m375_fixedOP[tmin:,anal_ind]+m375_fixedOP[tmin:,ss_ind]-m375_fixedOP[tmin:,sn_ind])/m375_fixedOP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[tmin:,x_ind], (m500_fixedSP[tmin:,DP_ind]-m500_fixedSP[tmin:,anal_ind]+m500_fixedSP[tmin:,ss_ind]-m500_fixedSP[tmin:,sn_ind])/m500_fixedSP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[tmin:,x_ind], (m500_bothfree[tmin:,DP_ind]-m500_bothfree[tmin:,anal_ind]+m500_bothfree[tmin:,ss_ind]-m500_bothfree[tmin:,sn_ind])/m500_bothfree[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[tmin:,x_ind], (m500_fixedOP[tmin:,DP_ind]-m500_fixedOP[tmin:,anal_ind]+m500_fixedOP[tmin:,ss_ind]-m500_fixedOP[tmin:,sn_ind])/m500_fixedOP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000plot
plt.scatter(m1000_fixedSP[tmin:,x_ind], (m1000_fixedSP[tmin:,DP_ind]-m1000_fixedSP[tmin:,anal_ind]+m1000_fixedSP[tmin:,ss_ind]-m1000_fixedSP[tmin:,sn_ind])/m1000_fixedSP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[tmin:,x_ind], (m1000_bothfree[tmin:,DP_ind]-m1000_bothfree[tmin:,anal_ind]+m1000_bothfree[tmin:,ss_ind]-m1000_bothfree[tmin:,sn_ind])/m1000_bothfree[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[tmin:,x_ind], (m1000_fixedOP[tmin:,DP_ind]-m1000_fixedOP[tmin:,anal_ind]+m1000_fixedOP[tmin:,ss_ind]-m1000_fixedOP[tmin:,sn_ind])/m1000_fixedOP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("full stress misfit/B",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-2,  2); 

# plot model DP misfit vs. x-axis variable

ax=fig.add_subplot(gs[0,1])

# 50
plt.scatter(m50_fixedSP[tmin:,x_ind], -1.0*(m50_fixedSP[tmin:,DP_ind]-m50_fixedSP[tmin:,anal_ind])/m50_fixedSP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[tmin:,x_ind], -1.0*(m50_bothfree[tmin:,DP_ind]-m50_bothfree[tmin:,anal_ind])/m50_bothfree[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[tmin:,x_ind], -1.0*(m50_fixedOP[tmin:,DP_ind]-m50_fixedOP[tmin:,anal_ind])/m50_fixedOP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[tmin:,x_ind], -1.0*(m250_fixedSP[tmin:,DP_ind]-m250_fixedSP[tmin:,anal_ind])/m250_fixedSP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[tmin:,x_ind], -1.0*(m250_bothfree[tmin:,DP_ind]-m250_bothfree[tmin:,anal_ind])/m250_bothfree[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[tmin:,x_ind], -1.0*(m250_fixedOP[tmin:,DP_ind]-m250_fixedOP[tmin:,anal_ind])/m250_fixedOP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[tmin:,x_ind], -1.0*(m375_fixedSP[tmin:,DP_ind]-m375_fixedSP[tmin:,anal_ind])/m375_fixedSP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[tmin:,x_ind], -1.0*(m375_bothfree[tmin:,DP_ind]-m375_bothfree[tmin:,anal_ind])/m375_bothfree[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[tmin:,x_ind], -1.0*(m375_fixedOP[tmin:,DP_ind]-m375_fixedOP[tmin:,anal_ind])/m375_fixedOP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[tmin:,x_ind], -1.0*(m500_fixedSP[tmin:,DP_ind]-m500_fixedSP[tmin:,anal_ind])/m500_fixedSP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[tmin:,x_ind], -1.0*(m500_bothfree[tmin:,DP_ind]-m500_bothfree[tmin:,anal_ind])/m500_bothfree[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[tmin:,x_ind], -1.0*(m500_fixedOP[tmin:,DP_ind]-m500_fixedOP[tmin:,anal_ind])/m500_fixedOP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[tmin:,x_ind], -1.0*(m1000_fixedSP[tmin:,DP_ind]-m1000_fixedSP[tmin:,anal_ind])/m1000_fixedSP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[tmin:,x_ind], -1.0*(m1000_bothfree[tmin:,DP_ind]-m1000_bothfree[tmin:,anal_ind])/m1000_bothfree[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[tmin:,x_ind], -1.0*(m1000_fixedOP[tmin:,DP_ind]-m1000_fixedOP[tmin:,anal_ind])/m1000_fixedOP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("-DP misfit [MPa]/B",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-2,  2); 


ax=fig.add_subplot(gs[1,0])

# 50
plt.scatter(m50_fixedSP[tmin:,x_ind], -1.0*m50_fixedSP[tmin:,sn_ind]/m50_fixedSP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[tmin:,x_ind],-1.0*m50_bothfree[tmin:,sn_ind]/m50_bothfree[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[tmin:,x_ind], -1.0*m50_fixedOP[tmin:,sn_ind]/m50_fixedOP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[tmin:,x_ind], -1.0*m250_fixedSP[tmin:,sn_ind]/m250_fixedSP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[tmin:,x_ind],-1.0*m250_bothfree[tmin:,sn_ind]/m250_bothfree[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[tmin:,x_ind], -1.0*m250_fixedOP[tmin:,sn_ind]/m250_fixedOP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[tmin:,x_ind], -1.0*m375_fixedSP[tmin:,sn_ind]/m375_fixedSP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[tmin:,x_ind],-1.0*m375_bothfree[tmin:,sn_ind]/m375_bothfree[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[tmin:,x_ind], -1.0*m375_fixedOP[tmin:,sn_ind]/m375_fixedOP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[tmin:,x_ind], -1.0*m500_fixedSP[tmin:,sn_ind]/m500_fixedSP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[tmin:,x_ind],-1.0*m500_bothfree[tmin:,sn_ind]/m500_bothfree[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[tmin:,x_ind], -1.0*m500_fixedOP[tmin:,sn_ind]/m500_fixedOP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[tmin:,x_ind], -1.0*m1000_fixedSP[tmin:,sn_ind]/m1000_fixedSP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[tmin:,x_ind],-1.0*m1000_bothfree[tmin:,sn_ind]/m1000_bothfree[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[tmin:,x_ind], -1.0*m1000_fixedOP[tmin:,sn_ind]/m1000_fixedOP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("-dQ/ds /B",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-2,  2); 

### plot normal stress term as a function of x-axis variable ###

ax=fig.add_subplot(gs[1,1])

# 50
plt.scatter(m50_fixedSP[tmin:,x_ind], m50_fixedSP[tmin:,ss_ind]/m50_fixedSP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m50_bothfree[tmin:,x_ind],m50_bothfree[tmin:,ss_ind]/m50_bothfree[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m50_fixedOP[tmin:,x_ind], m50_fixedOP[tmin:,ss_ind]/m50_fixedOP[tmin:,anal_ind],s=10,color='tan',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 250
plt.scatter(m250_fixedSP[tmin:,x_ind], m250_fixedSP[tmin:,ss_ind]/m250_fixedSP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m250_bothfree[tmin:,x_ind],m250_bothfree[tmin:,ss_ind]/m250_bothfree[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m250_fixedOP[tmin:,x_ind], m250_fixedOP[tmin:,ss_ind]/m250_fixedOP[tmin:,anal_ind],s=10,color='peru',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 375
plt.scatter(m375_fixedSP[tmin:,x_ind], m375_fixedSP[tmin:,ss_ind]/m375_fixedSP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m375_bothfree[tmin:,x_ind],m375_bothfree[tmin:,ss_ind]/m375_bothfree[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m375_fixedOP[tmin:,x_ind], m375_fixedOP[tmin:,ss_ind]/m375_fixedOP[tmin:,anal_ind],s=10,color='firebrick',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 500
plt.scatter(m500_fixedSP[tmin:,x_ind], m500_fixedSP[tmin:,ss_ind]/m500_fixedSP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m500_bothfree[tmin:,x_ind],m500_bothfree[tmin:,ss_ind]/m500_bothfree[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m500_fixedOP[tmin:,x_ind], m500_fixedOP[tmin:,ss_ind]/m500_fixedOP[tmin:,anal_ind],s=10,color='maroon',edgecolor='black',linewidth=0.25,zorder=3,marker='^')
# 1000
plt.scatter(m1000_fixedSP[tmin:,x_ind], m1000_fixedSP[tmin:,ss_ind]/m1000_fixedSP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='v')
plt.scatter(m1000_bothfree[tmin:,x_ind],m1000_bothfree[tmin:,ss_ind]/m1000_bothfree[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3)
plt.scatter(m1000_fixedOP[tmin:,x_ind], m1000_fixedOP[tmin:,ss_ind]/m1000_fixedOP[tmin:,anal_ind],s=10,color='black',edgecolor='black',linewidth=0.25,zorder=3,marker='^')

plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel("KN/B",size=6)
plt.xlabel(x_label,size=6)
plt.ylim(-2,  2); 


plt.subplots_adjust(wspace=0.5)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
# plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')
