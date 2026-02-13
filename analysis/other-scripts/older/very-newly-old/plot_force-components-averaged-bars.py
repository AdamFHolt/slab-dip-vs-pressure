#!/bin/python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions import get_avg_forces_nondim, get_avg_forces_nondim_curvethresh
import matplotlib.font_manager as fm

import matplotlib.font_manager as fm
font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
myriad_pro = fm.FontProperties(fname=font_path)

# Now, set it globally
mpl.rcParams['font.family'] = 'Myriad Pro'  # Now it should work if properly installed!
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['axes.labelpad'] = 1.25
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['xtick.major.pad'] = 2
mpl.rcParams['ytick.major.pad'] = 2
mpl.rcParams['xtick.major.size'] = 3.5
mpl.rcParams['ytick.major.size'] = 3.5
mpl.rcParams['xtick.minor.size'] = 1.75
mpl.rcParams['ytick.minor.size'] = 1.75

analysis_depth  = float(sys.argv[1]) 
analysis_depth_dz = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

tactual_min = 11 # first time step to use
tmin = tactual_min - 8

plot_name_png = ''.join(['plots/DP-comparisons/compilations/force-components-averaged-bars.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/force-components-averaged-bars.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.pdf'])

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


# col 1 = visc, col 2 = DP anal, col 3 = full mod, col 4 = DP, col 5 = shear, col 6 = norm, col 6 = K
fixedSP_array = np.zeros((5,9))
fixedSP_array[0:5,0] = [50,250,375,500,700]
curve_thresh = 0.0015 # [1/km]

fixedSP_array[0,1], fixedSP_array[0,2], fixedSP_array[0,3], fixedSP_array[0,4], fixedSP_array[0,5], fixedSP_array[0,6], fixedSP_array[0,7], fixedSP_array[0,8] = get_avg_forces_nondim_curvethresh(m50_fixedSP,tmin,curve_thresh)
fixedSP_array[1,1], fixedSP_array[1,2], fixedSP_array[1,3], fixedSP_array[1,4], fixedSP_array[1,5], fixedSP_array[1,6], fixedSP_array[1,7], fixedSP_array[1,8]  = get_avg_forces_nondim_curvethresh(m250_fixedSP,tmin,curve_thresh)
fixedSP_array[2,1], fixedSP_array[2,2], fixedSP_array[2,3], fixedSP_array[2,4], fixedSP_array[2,5], fixedSP_array[2,6], fixedSP_array[2,7], fixedSP_array[2,8] = get_avg_forces_nondim_curvethresh(m375_fixedSP,tmin,curve_thresh)
fixedSP_array[3,1], fixedSP_array[3,2], fixedSP_array[3,3], fixedSP_array[3,4], fixedSP_array[3,5], fixedSP_array[3,6], fixedSP_array[3,7], fixedSP_array[3,8] = get_avg_forces_nondim_curvethresh(m500_fixedSP,tmin,curve_thresh)
fixedSP_array[4,1], fixedSP_array[4,2], fixedSP_array[4,3], fixedSP_array[4,4], fixedSP_array[4,5], fixedSP_array[4,6], fixedSP_array[4,7], fixedSP_array[4,8] = get_avg_forces_nondim_curvethresh(m1000_fixedSP,tmin,curve_thresh)

bothfree_array = np.zeros((5,9))
bothfree_array[0:5,0] = [50,250,375,500,700]
bothfree_array[0,1], bothfree_array[0,2], bothfree_array[0,3], bothfree_array[0,4], bothfree_array[0,5], bothfree_array[0,6], bothfree_array[0,7], bothfree_array[0,8] = get_avg_forces_nondim_curvethresh(m50_bothfree,tmin,curve_thresh)
bothfree_array[1,1], bothfree_array[1,2], bothfree_array[1,3], bothfree_array[1,4], bothfree_array[1,5], bothfree_array[1,6], bothfree_array[1,7], bothfree_array[1,8] = get_avg_forces_nondim_curvethresh(m250_bothfree,tmin,curve_thresh)
bothfree_array[2,1], bothfree_array[2,2], bothfree_array[2,3], bothfree_array[2,4], bothfree_array[2,5], bothfree_array[2,6], bothfree_array[2,7], bothfree_array[2,8] = get_avg_forces_nondim_curvethresh(m375_bothfree,tmin,curve_thresh)
bothfree_array[3,1], bothfree_array[3,2], bothfree_array[3,3], bothfree_array[3,4], bothfree_array[3,5], bothfree_array[3,6], bothfree_array[3,7], bothfree_array[3,8] = get_avg_forces_nondim_curvethresh(m500_bothfree,tmin,curve_thresh)
bothfree_array[4,1], bothfree_array[4,2], bothfree_array[4,3], bothfree_array[4,4], bothfree_array[4,5], bothfree_array[4,6], bothfree_array[4,7], bothfree_array[4,8] = get_avg_forces_nondim_curvethresh(m1000_bothfree,tmin,curve_thresh)

fixedOP_array = np.zeros((5,9))
fixedOP_array[0:5,0] = [50,250,375,500,700]
fixedOP_array[0,1], fixedOP_array[0,2], fixedOP_array[0,3], fixedOP_array[0,4], fixedOP_array[0,5], fixedOP_array[0,6], fixedOP_array[0,7], fixedOP_array[0,8]  = get_avg_forces_nondim_curvethresh(m50_fixedOP,tmin,curve_thresh)
fixedOP_array[1,1], fixedOP_array[1,2], fixedOP_array[1,3], fixedOP_array[1,4], fixedOP_array[1,5], fixedOP_array[1,6], fixedOP_array[1,7], fixedOP_array[1,8]  = get_avg_forces_nondim_curvethresh(m250_fixedOP,tmin,curve_thresh)
fixedOP_array[2,1], fixedOP_array[2,2], fixedOP_array[2,3], fixedOP_array[2,4], fixedOP_array[2,5], fixedOP_array[2,6], fixedOP_array[2,7], fixedOP_array[2,8]  = get_avg_forces_nondim_curvethresh(m375_fixedOP,tmin,curve_thresh)
fixedOP_array[3,1], fixedOP_array[3,2], fixedOP_array[3,3], fixedOP_array[3,4], fixedOP_array[3,5], fixedOP_array[3,6], fixedOP_array[3,7], fixedOP_array[3,8]  = get_avg_forces_nondim_curvethresh(m500_fixedOP,tmin,curve_thresh)
fixedOP_array[4,1], fixedOP_array[4,2], fixedOP_array[4,3], fixedOP_array[4,4], fixedOP_array[4,5], fixedOP_array[4,6], fixedOP_array[4,7], fixedOP_array[4,8]  = get_avg_forces_nondim_curvethresh(m1000_fixedOP,tmin,curve_thresh)

fig=plt.figure()
gs=GridSpec(10,5) 

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


ax=fig.add_subplot(gs[1:5,0:4])

shift=35
w=25

# fixed SP models
for i in range(len(fixedSP_array)):

    # start with in-slab stresses 
    if (fixedSP_array[i,5]/fixedSP_array[i,1]) < 0 and (fixedSP_array[i,4]/fixedSP_array[i,1]) < 0:
        bottom1  = (fixedSP_array[i,4]+fixedSP_array[i,5])/fixedSP_array[i,1]
        thick1   = -1.0*fixedSP_array[i,4]/fixedSP_array[i,1]
        bottom2  =  (fixedSP_array[i,5])/fixedSP_array[i,1]
        thick2   = -1.0*fixedSP_array[i,5]/fixedSP_array[i,1]
        #
        plt.bar(fixedSP_array[i,0]-shift, thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(fixedSP_array[i,0]-shift, thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    elif (fixedSP_array[i,5]/fixedSP_array[i,1]) >= 0 and (fixedSP_array[i,4]/fixedSP_array[i,1]) < 0:
        bottom1  = fixedSP_array[i,4]/fixedSP_array[i,1]
        thick1   = -1.0*fixedSP_array[i,4]/fixedSP_array[i,1]
        bottom2  =  0
        thick2   = fixedSP_array[i,5]/fixedSP_array[i,1]
        #
        plt.bar(fixedSP_array[i,0]-shift, thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(fixedSP_array[i,0]-shift, thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    elif (fixedSP_array[i,5]/fixedSP_array[i,1]) < 0 and (fixedSP_array[i,4]/fixedSP_array[i,1]) >= 0:
        bottom1  = fixedSP_array[i,5]/fixedSP_array[i,1]
        thick1   = -1.0*fixedSP_array[i,5]/fixedSP_array[i,1]
        bottom2  =  0
        thick2   = fixedSP_array[i,4]/fixedSP_array[i,1]
        #
        plt.bar(fixedSP_array[i,0]-shift, thick1, bottom=bottom1, width=w, color='cornflowerblue', zorder=3)
        plt.bar(fixedSP_array[i,0]-shift, thick2, bottom=bottom2, width=w, color='mediumblue', zorder=3)

    elif (fixedSP_array[i,5]/fixedSP_array[i,1]) and (fixedSP_array[i,4]/fixedSP_array[i,1]) > 0:
        bottom1  = 0
        thick1   = fixedSP_array[i,4]/fixedSP_array[i,1]
        bottom2  =  fixedSP_array[i,4]/fixedSP_array[i,1]
        thick2   = fixedSP_array[i,5]/fixedSP_array[i,1]
        #
        plt.bar(fixedSP_array[i,0]-shift, thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(fixedSP_array[i,0]-shift, thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)
    

    # now DP
    if fixedSP_array[i,3]/fixedSP_array[i,1] > 0:
        bottom3  =  bottom2 + thick2
        thick3   = fixedSP_array[i,3]/fixedSP_array[i,1]
        plt.bar(fixedSP_array[i,0]-shift, thick3, bottom=bottom3, width=w, color='red', zorder=3)
    else:
        bottom3  =  (bottom2 + thick2) + (fixedSP_array[i,3]/fixedSP_array[i,1])
        thick3   = -1.0*(fixedSP_array[i,3]/fixedSP_array[i,1])
        plt.bar(fixedSP_array[i,0]-shift, thick3, bottom=bottom3, width=0.8*w, color='red', zorder=3)


# total model force
plt.scatter(fixedSP_array[:,0]-shift,fixedSP_array[:,2]/fixedSP_array[:,1],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='v')
# stdev of total model force
plt.plot([fixedSP_array[:,0]-shift,fixedSP_array[:,0]-shift], [(fixedSP_array[:,2]/fixedSP_array[:,1])-fixedSP_array[:,6],(fixedSP_array[:,2]/fixedSP_array[:,1])+fixedSP_array[:,6]],  color='silver',linewidth=1,zorder=5)


# free SP models
for i in range(len(bothfree_array)):

    # start with in-slab stresses 
    if (bothfree_array[i,5]/bothfree_array[i,1]) < 0 and (bothfree_array[i,4]/bothfree_array[i,1]) < 0:
        bottom1  = (bothfree_array[i,4]+bothfree_array[i,5])/bothfree_array[i,1]
        thick1   = -1.0*bothfree_array[i,4]/bothfree_array[i,1]
        bottom2  =  (bothfree_array[i,5])/bothfree_array[i,1]
        thick2   = -1.0*bothfree_array[i,5]/bothfree_array[i,1]
        #
        plt.bar(bothfree_array[i,0], thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(bothfree_array[i,0], thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    elif (bothfree_array[i,5]/bothfree_array[i,1]) >= 0 and (bothfree_array[i,4]/bothfree_array[i,1]) < 0:
        bottom1  = bothfree_array[i,4]/bothfree_array[i,1]
        thick1   = -1.0*bothfree_array[i,4]/bothfree_array[i,1]
        bottom2  =  0
        thick2   = bothfree_array[i,5]/bothfree_array[i,1]
        #
        plt.bar(bothfree_array[i,0], thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(bothfree_array[i,0], thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    elif (bothfree_array[i,5]/bothfree_array[i,1]) < 0 and (bothfree_array[i,4]/bothfree_array[i,1]) >= 0:
        bottom1  = bothfree_array[i,5]/bothfree_array[i,1]
        thick1   = -1.0*bothfree_array[i,5]/bothfree_array[i,1]
        bottom2  =  0
        thick2   = bothfree_array[i,4]/bothfree_array[i,1]
        #
        plt.bar(bothfree_array[i,0], thick1, bottom=bottom1, width=w, color='cornflowerblue', zorder=3)
        plt.bar(bothfree_array[i,0], thick2, bottom=bottom2, width=w, color='mediumblue', zorder=3)

    elif (bothfree_array[i,5]/bothfree_array[i,1]) and (bothfree_array[i,4]/bothfree_array[i,1]) > 0:
        bottom1  = 0
        thick1   = bothfree_array[i,4]/bothfree_array[i,1]
        bottom2  =  bothfree_array[i,4]/bothfree_array[i,1]
        thick2   = bothfree_array[i,5]/bothfree_array[i,1]
        #
        plt.bar(bothfree_array[i,0], thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(bothfree_array[i,0], thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    # now DP
    if bothfree_array[i,3]/bothfree_array[i,1] > 0:
        bottom3  =  bottom2 + thick2
        thick3   = bothfree_array[i,3]/bothfree_array[i,1]
        plt.bar(bothfree_array[i,0], thick3, bottom=bottom3, width=w, color='red', zorder=3)
    else:
        bottom3  =  (bottom2 + thick2) + (bothfree_array[i,3]/bothfree_array[i,1])
        thick3   = -1.0*(bothfree_array[i,3]/bothfree_array[i,1])
        plt.bar(bothfree_array[i,0], thick3, bottom=bottom3, width=0.8*w, color='red', zorder=3)

# total model force
plt.scatter(bothfree_array[:,0],bothfree_array[:,2]/bothfree_array[:,1],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')
# stdev of total model force
plt.plot([bothfree_array[:,0],bothfree_array[:,0]], [(bothfree_array[:,2]/bothfree_array[:,1])-bothfree_array[:,6],(bothfree_array[:,2]/bothfree_array[:,1])+bothfree_array[:,6]],  color='silver',linewidth=1,zorder=5)


# fixed OP models
for i in range(len(fixedOP_array)):

    if fixedOP_array[i,4] == 0:
        break

    # start with in-slab stresses 
    if (fixedOP_array[i,5]/fixedOP_array[i,1]) < 0 and (fixedOP_array[i,4]/fixedOP_array[i,1]) < 0:
        bottom1  = (fixedOP_array[i,4]+fixedOP_array[i,5])/fixedOP_array[i,1]
        thick1   = -1.0*fixedOP_array[i,4]/fixedOP_array[i,1]
        bottom2  =  (fixedOP_array[i,5])/fixedOP_array[i,1]
        thick2   = -1.0*fixedOP_array[i,5]/fixedOP_array[i,1]
        #
        plt.bar(fixedOP_array[i,0]+shift, thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(fixedOP_array[i,0]+shift, thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    elif (fixedOP_array[i,5]/fixedOP_array[i,1]) >= 0 and (fixedOP_array[i,4]/fixedOP_array[i,1]) < 0:
        bottom1  = fixedOP_array[i,4]/fixedOP_array[i,1]
        thick1   = -1.0*fixedOP_array[i,4]/fixedOP_array[i,1]
        bottom2  =  0
        thick2   = fixedOP_array[i,5]/fixedOP_array[i,1]
        #
        plt.bar(fixedOP_array[i,0]+shift, thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(fixedOP_array[i,0]+shift, thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    elif (fixedOP_array[i,5]/fixedOP_array[i,1]) < 0 and (fixedOP_array[i,4]/fixedOP_array[i,1]) >= 0:
        bottom1  = fixedOP_array[i,5]/fixedOP_array[i,1]
        thick1   = -1.0*fixedOP_array[i,5]/fixedOP_array[i,1]
        bottom2  =  0
        thick2   = fixedOP_array[i,4]/fixedOP_array[i,1]
        #
        plt.bar(fixedOP_array[i,0]+shift, thick1, bottom=bottom1, width=w, color='cornflowerblue', zorder=3)
        plt.bar(fixedOP_array[i,0]+shift, thick2, bottom=bottom2, width=w, color='mediumblue', zorder=3)

    elif (fixedOP_array[i,5]/fixedOP_array[i,1]) and (fixedOP_array[i,4]/fixedOP_array[i,1]) > 0:
        bottom1  = 0
        thick1   = fixedOP_array[i,4]/fixedOP_array[i,1]
        bottom2  =  fixedOP_array[i,4]/fixedOP_array[i,1]
        thick2   = fixedOP_array[i,5]/fixedOP_array[i,1]
        #
        plt.bar(fixedOP_array[i,0]+shift, thick1, bottom=bottom1, width=w, color='mediumblue', zorder=3)
        plt.bar(fixedOP_array[i,0]+shift, thick2, bottom=bottom2, width=w, color='cornflowerblue', zorder=3)

    # now DP
    if fixedOP_array[i,3]/fixedOP_array[i,1] > 0:
        bottom3  =  bottom2 + thick2
        thick3   = fixedOP_array[i,3]/fixedOP_array[i,1]
        plt.bar(fixedOP_array[i,0]+shift, thick3, bottom=bottom3, width=w, color='red', zorder=3)
    else:
        bottom3  =  (bottom2 + thick2) + (fixedOP_array[i,3]/fixedOP_array[i,1])
        thick3   = -1.0*(fixedOP_array[i,3]/fixedOP_array[i,1])
        plt.bar(fixedOP_array[i,0]+shift, thick3, bottom=bottom3, width=0.8*w, color='red', zorder=3)

# total model force
for i in range(len(fixedOP_array)):
    if fixedOP_array[i,4] == 0:
        pass
    else:
        plt.scatter(fixedOP_array[i,0]+shift,fixedOP_array[i,2]/fixedOP_array[i,1],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='^')
# stdev of total model force
for i in range(len(fixedOP_array)):
    if fixedOP_array[i,4] == 0:
        pass
    else:
        plt.plot([fixedOP_array[i,0]+shift,fixedOP_array[i,0]+shift], [(fixedOP_array[i,2]/fixedOP_array[i,1])-fixedOP_array[i,6],(fixedOP_array[i,2]/fixedOP_array[i,1])+fixedOP_array[i,6]],  color='silver',linewidth=1,zorder=5)

ax.set_xticks( [50,250, 375, 500,725] )
ax.set_xticklabels( ['50','250','375','500','1000'] )
plt.axhline(y=0, color='silver',linestyle='--',linewidth=1, zorder=1)
plt.axhline(y=1, color='gray',linestyle='-',linewidth=1.25, zorder=1)
plt.ylabel(r'$F$ / $B_{slab}$')
plt.xlabel(r'$\eta_{slab}$ / $\eta_{mantle}$')
plt.ylim(-0.5,  1.2);
plt.xlim(-50,825);
#ax.set_yticks([-0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
#ax.set_yticklabels([-0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
ax.set_yticks([0,0.5,1.0])
ax.set_yticklabels([0,0.5,1.0])
#ax.set_yticklabels([-0.4,-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

plt.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)


ax2=fig.add_subplot(gs[0,0:4])

# average curvature
#plt.scatter(fixedSP_array[:,0]-shift,fixedSP_array[:,8],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=5,marker='v')
#plt.scatter(bothfree_array[:,0],bothfree_array[:,8],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=5,marker='o')
#plt.scatter(fixedOP_array[:,0]+shift,fixedOP_array[:,8],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=5,marker='^')
plt.scatter(fixedSP_array[:,0]-shift,fixedSP_array[:,7],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='v')
plt.scatter(bothfree_array[:,0],bothfree_array[:,7],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='o')
for i in range(len(fixedOP_array)):
    if fixedOP_array[i,7] == 0:
        pass
    else:
        plt.scatter(fixedOP_array[i,0]+shift,fixedOP_array[i,7],s=30,color='black',edgecolor='black',linewidth=0.25,zorder=6,marker='^')

ax2.axhline(y=0, color='silver',linestyle='--',linewidth=1, zorder=0)
ax2.set_xticks( [50,250, 375, 500,725] )
ax2.set_xticklabels( [] )
plt.ylabel(r'$K$  [1/km]')
plt.ylim(-0.0005,0.0025);
plt.xlim(-50,825);


ax2.minorticks_on()
ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.0005))



plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')



