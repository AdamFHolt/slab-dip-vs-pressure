#!/bin/python
import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions import get_misfit_mean_and_stdev, get_curvature_mean_and_stdev, get_misfit_mean_and_stdev_nondim
from functions_plotting import plot_forcecomponent_fullstressmisfit, plot_forcecomponent_dqds, plot_forcecomponent_dpmisfit
from functions_plotting import plot_forcecomponent_fullstressmisfit_overturned, plot_forcecomponent_dqds_overturned, plot_forcecomponent_dpmisfit_overturned
import matplotlib.font_manager as fm
font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
myriad_pro = fm.FontProperties(fname=font_path)

mpl.rcParams['font.family'] = 'Myriad Pro'  # Now it should work if properly installed!
mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['axes.labelpad'] = 1.25
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['xtick.major.pad'] = 2
mpl.rcParams['ytick.major.pad'] = 2
mpl.rcParams['xtick.major.size'] = 2.5
mpl.rcParams['ytick.major.size'] = 2.5
mpl.rcParams['xtick.minor.size'] = 1.25
mpl.rcParams['ytick.minor.size'] = 1.25

analysis_depth  = float(sys.argv[1]) 
analysis_depth_dz = float(sys.argv[2])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[3])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[4])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)
coeff = 100./1624.0


tactual_min = 11 # first time step to use
tmin = tactual_min - 8

plot_name_png = ''.join(['plots/DP-comparisons/compilations/forces-vs-scaling.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.no-ot.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/forces-vs-scaling.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.no-ot.pdf'])

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
text50_bothfree 		= ''.join(['text_files/TESTB/',name1_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedSP  		= ''.join(['text_files/TESTB/',name1_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedOP  		= ''.join(['text_files/TESTB/',name1_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_bothfree    	= ''.join(['text_files/TESTB/',name3_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedSP     	= ''.join(['text_files/TESTB/',name3_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedOP		    = ''.join(['text_files/TESTB/',name3_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_bothfree 		= ''.join(['text_files/TESTB/',name4_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedSP  		= ''.join(['text_files/TESTB/',name4_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedOP  		= ''.join(['text_files/TESTB/',name4_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_bothfree		= ''.join(['text_files/TESTB/',name5_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedSP  		= ''.join(['text_files/TESTB/',name5_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedOP  		= ''.join(['text_files/TESTB/',name5_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_bothfree    	= ''.join(['text_files/TESTB/',name7_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedSP     	= ''.join(['text_files/TESTB/',name7_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedOP	    	= ''.join(['text_files/TESTB/',name7_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

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
vs_ind   = 20

x_ind = vc_ind
x_label = "V K visc [Pa]"

s_to_yr = 1./(365.25*24*3600)
cmyr_to_ms = 0.01/(365.25*24*3600)

mant_visc = 2.5e20
m50_bothfree[:,vc_ind] 	=  coeff * m50_bothfree[:,K_ind]    * 50 * mant_visc   * m50_bothfree[:,vc_ind]   * cmyr_to_ms * 1e-6
m50_fixedSP[:,vc_ind] 	=  coeff * m50_fixedSP[:,K_ind]     * 50 * mant_visc   * m50_fixedSP[:,vc_ind]    * cmyr_to_ms * 1e-6
m50_fixedOP[:,vc_ind] 	=  coeff * m50_fixedOP[:,K_ind]     * 50 * mant_visc   * m50_fixedOP[:,vc_ind]    * cmyr_to_ms * 1e-6
m250_bothfree[:,vc_ind] 	=  coeff * m250_bothfree[:,K_ind]   * 250 * mant_visc  * m250_bothfree[:,vc_ind]  * cmyr_to_ms * 1e-6
m250_fixedSP[:,vc_ind] 	=  coeff * m250_fixedSP[:,K_ind]    * 250 * mant_visc  * m250_fixedSP[:,vc_ind]   * cmyr_to_ms * 1e-6
m250_fixedOP[:,vc_ind] 	=  coeff * m250_fixedOP[:,K_ind]    * 250 * mant_visc  * m250_fixedOP[:,vc_ind]   * cmyr_to_ms * 1e-6
m375_bothfree[:,vc_ind] 	=  coeff * m375_bothfree[:,K_ind]   * 375 * mant_visc  * m375_bothfree[:,vc_ind]  * cmyr_to_ms * 1e-6
m375_fixedSP[:,vc_ind] 	=  coeff * m375_fixedSP[:,K_ind]    * 375 * mant_visc  * m375_fixedSP[:,vc_ind]   * cmyr_to_ms * 1e-6
m375_fixedOP[:,vc_ind] 	=  coeff * m375_fixedOP[:,K_ind]    * 375 * mant_visc  * m375_fixedOP[:,vc_ind]   * cmyr_to_ms * 1e-6
m500_bothfree[:,vc_ind] 	=  coeff * m500_bothfree[:,K_ind]   * 500 * mant_visc  * m500_bothfree[:,vc_ind]  * cmyr_to_ms * 1e-6
m500_fixedSP[:,vc_ind] 	=  coeff * m500_fixedSP[:,K_ind]    * 500 * mant_visc  * m500_fixedSP[:,vc_ind]   * cmyr_to_ms * 1e-6
m500_fixedOP[:,vc_ind] 	=  coeff * m500_fixedOP[:,K_ind]    * 500 * mant_visc  * m500_fixedOP[:,vc_ind]   * cmyr_to_ms * 1e-6
m1000_bothfree[:,vc_ind] =  coeff * m1000_bothfree[:,K_ind]  * 1000 * mant_visc * m1000_bothfree[:,vc_ind] * cmyr_to_ms * 1e-6
m1000_fixedSP[:,vc_ind] 	=  coeff * m1000_fixedSP[:,K_ind]   * 1000 * mant_visc * m1000_fixedSP[:,vc_ind]  * cmyr_to_ms * 1e-6
m1000_fixedOP[:,vc_ind] 	=  coeff * m1000_fixedOP[:,K_ind]   * 1000 *mant_visc  *  m1000_fixedOP[:,vc_ind] * cmyr_to_ms * 1e-6

gs=GridSpec(1,3)

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

curve_thresh = 1e9
misfit_color = 'gold'
plot_forcecomponent_fullstressmisfit(tmin,m50_fixedSP,curve_thresh,'tan','v',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m50_bothfree,curve_thresh,'tan','o',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m50_fixedOP,curve_thresh,'tan','^',misfit_color)

plot_forcecomponent_fullstressmisfit(tmin,m250_fixedSP,curve_thresh,'peru','v',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m250_bothfree,curve_thresh,'peru','o',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m250_fixedOP,curve_thresh,'peru','^',misfit_color)

plot_forcecomponent_fullstressmisfit(tmin,m375_fixedSP,curve_thresh,'firebrick','v',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m375_bothfree,curve_thresh,'firebrick','o',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m375_fixedOP,curve_thresh,'firebrick','^',misfit_color)

plot_forcecomponent_fullstressmisfit(tmin,m500_fixedSP,curve_thresh,'maroon','v',misfit_color)
plot_forcecomponent_fullstressmisfit(tmin,m500_bothfree,curve_thresh,'maroon','o',misfit_color)
plot_forcecomponent_fullstressmisfit_overturned(tmin,m500_fixedOP,'maroon','^')

plot_forcecomponent_fullstressmisfit(tmin,m1000_fixedSP,curve_thresh,'black','v',misfit_color)
plot_forcecomponent_fullstressmisfit_overturned(tmin,m1000_bothfree,'black','o')
plot_forcecomponent_fullstressmisfit_overturned(tmin,m1000_fixedOP,'black','^')

# axis stuff
plt.ylim(-10,  17.5); 
plt.xlim(-0.001,  0.003); 
plt.ylabel("full stress misfit [MPa]")
plt.ylabel(r'$(\Delta P  +  \sigma_{slab}) - B_{slab}$   [MPa]')
plt.xlabel(r'$K$    [1/km]')
ax.set_xticks([-0.001, 0, 0.001, 0.002, 0.003])
ax.set_xticklabels([-0.001, 0, 0.001, 0.002, 0.003])
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.0005))
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
plt.grid(True, which='both', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)


# plot model DP misfit vs. x-axis variable

ax=fig.add_subplot(gs[0,2])

plot_forcecomponent_dpmisfit(tmin,m50_fixedSP,curve_thresh,x_ind,'tan','v',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m50_bothfree,curve_thresh,x_ind,'tan','o',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m50_fixedOP,curve_thresh,x_ind,'tan','^',misfit_color)

plot_forcecomponent_dpmisfit(tmin,m250_fixedSP,curve_thresh,x_ind,'peru','v',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m250_bothfree,curve_thresh,x_ind,'peru','o',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m250_fixedOP,curve_thresh,x_ind,'peru','^',misfit_color)

plot_forcecomponent_dpmisfit(tmin,m375_fixedSP,curve_thresh,x_ind,'firebrick','v',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m375_bothfree,curve_thresh,x_ind,'firebrick','o',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m375_fixedOP,curve_thresh,x_ind,'firebrick','^',misfit_color)

plot_forcecomponent_dpmisfit(tmin,m500_fixedSP,curve_thresh,x_ind,'maroon','v',misfit_color)
plot_forcecomponent_dpmisfit(tmin,m500_bothfree,curve_thresh,x_ind,'maroon','o',misfit_color)
plot_forcecomponent_dpmisfit_overturned(tmin,m500_fixedOP,x_ind,'maroon','^')

plot_forcecomponent_dpmisfit(tmin,m1000_fixedSP,curve_thresh,x_ind,'black','v',misfit_color)
plot_forcecomponent_dpmisfit_overturned(tmin,m1000_bothfree,x_ind,'black','o')
plot_forcecomponent_dpmisfit_overturned(tmin,m1000_fixedOP,x_ind,'black','^')

# axis stuff
plt.ylim(-10,  17.5);
plt.xlim(-5, 40);
plt.ylabel(r'$-(\Delta P - B_{slab}$)   [MPa]')
plt.xlabel(r'($\eta H K V_{C}$)/$L_{eff}$   [MPa]')
ax.set_xticks([0, 10, 20, 30, 40])
ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
plt.grid(True, which='both', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
x_values = np.linspace(-5, 40, 100)
plt.plot(x_values, x_values, color='gray', linestyle='--', linewidth=1, zorder=-1)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)

# plot dQ/dS vs. x-axis variable
ax=fig.add_subplot(gs[0,1])

plot_forcecomponent_dqds(tmin,m50_fixedSP,curve_thresh,x_ind,'tan','v',misfit_color)
plot_forcecomponent_dqds(tmin,m50_bothfree,curve_thresh,x_ind,'tan','o',misfit_color)
plot_forcecomponent_dqds(tmin,m50_fixedOP,curve_thresh,x_ind,'tan','^',misfit_color)

plot_forcecomponent_dqds(tmin,m250_fixedSP,curve_thresh,x_ind,'peru','v',misfit_color)
plot_forcecomponent_dqds(tmin,m250_bothfree,curve_thresh,x_ind,'peru','o',misfit_color)
plot_forcecomponent_dqds(tmin,m250_fixedOP,curve_thresh,x_ind,'peru','^',misfit_color)

plot_forcecomponent_dqds(tmin,m375_fixedSP,curve_thresh,x_ind,'firebrick','v',misfit_color)
plot_forcecomponent_dqds(tmin,m375_bothfree,curve_thresh,x_ind,'firebrick','o',misfit_color)
plot_forcecomponent_dqds(tmin,m375_fixedOP,curve_thresh,x_ind,'firebrick','^',misfit_color)

plot_forcecomponent_dqds(tmin,m500_fixedSP,curve_thresh,x_ind,'maroon','v',misfit_color)
plot_forcecomponent_dqds(tmin,m500_bothfree,curve_thresh,x_ind,'maroon','o',misfit_color)
plot_forcecomponent_dqds_overturned(tmin,m500_fixedOP,x_ind,'maroon','^')

plot_forcecomponent_dqds(tmin,m1000_fixedSP,curve_thresh,x_ind,'black','v',misfit_color)
plot_forcecomponent_dqds_overturned(tmin,m1000_bothfree,x_ind,'black','o')
plot_forcecomponent_dqds_overturned(tmin,m1000_fixedOP,x_ind,'black','^')

# axis stuff
# REGULAR
plt.ylim(-10,  17.5);
plt.xlim(-5, 40);
ax.set_xticks([0, 10, 20, 30, 40])

plt.ylabel(r'$-\frac{dQ}{ds}$   [MPa]')
plt.xlabel(r'($\eta H K V_{C}$)/$L_{eff}$   [MPa]')
ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
plt.grid(True, which='both', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
# Add a 1:1 line in the background
x_values = np.linspace(-5, 40, 100)
plt.plot(x_values, x_values, color='gray', linestyle='--', linewidth=1, zorder=-1)
plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)


# ### plot normal stress term as a function of x-axis variable ###
# ax=fig.add_subplot(gs[1,1])


# plot_forcecomponent_KN(tmin,m50_fixedSP,curve_thresh,'tan','v',misfit_color)
# plot_forcecomponent_KN(tmin,m50_bothfree,curve_thresh,'tan','o',misfit_color)
# plot_forcecomponent_KN(tmin,m50_fixedOP,curve_thresh,'tan','^',misfit_color)

# plot_forcecomponent_KN(tmin,m250_fixedSP,curve_thresh,'peru','v',misfit_color)
# plot_forcecomponent_KN(tmin,m250_bothfree,curve_thresh,'peru','o',misfit_color)
# plot_forcecomponent_KN(tmin,m250_fixedOP,curve_thresh,'peru','^',misfit_color)

# plot_forcecomponent_KN(tmin,m375_fixedSP,curve_thresh,'firebrick','v',misfit_color)
# plot_forcecomponent_KN(tmin,m375_bothfree,curve_thresh,'firebrick','o',misfit_color)
# plot_forcecomponent_KN(tmin,m375_fixedOP,curve_thresh,'firebrick','^',misfit_color)

# plot_forcecomponent_KN(tmin,m500_fixedSP,curve_thresh,'maroon','v',misfit_color)
# plot_forcecomponent_KN(tmin,m500_bothfree,curve_thresh,'maroon','o',misfit_color)
# plot_forcecomponent_KN_overturned(tmin,m500_fixedOP,curve_thresh,'maroon','^',misfit_color)

# plot_forcecomponent_KN(tmin,m1000_fixedSP,curve_thresh,'black','v',misfit_color)
# plot_forcecomponent_KN_overturned(tmin,m1000_bothfree,curve_thresh,'black','o',misfit_color)
# plot_forcecomponent_KN_overturned(tmin,m1000_fixedOP,curve_thresh,'black','^',misfit_color)


# # axis stuff
# plt.ylim(-3,  3); 
# plt.xlim(-0.001,  0.003); 
# plt.ylabel(r"$KN$   [MPa]")
# plt.xlabel(r'$K$    [1/km]')
# ax.set_xticks([-0.001, 0, 0.001, 0.002, 0.003])
# ax.set_xticklabels([-0.001, 0, 0.001, 0.002, 0.003])
# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.0005))
# ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
# plt.grid(True, which='both', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
# plt.axhline(y=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
# plt.axvline(x=0, color='lightgray',linestyle='-',linewidth=1, zorder=0)
# fixed_aspect_ratio(1)

plt.subplots_adjust(wspace=0.5)

plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=600)
plt.savefig(plot_name_pdf, format='pdf', bbox_inches='tight')
