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
from functions_plotting import plot_forcecomponent_fullstressmisfit, plot_forcecomponent_KN 
from functions_plotting import plot_forcecomponent_dqds, plot_forcecomponent_dpmisfit
from functions_plotting import append_DPcomponents
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

tactual_min = 11 # first time step to use
tmin = tactual_min - 8

plot_name_png = ''.join(['plots/DP-comparisons/compilations/Psp-DP-ratio.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/Psp-DP-ratio.z',str(analysis_depth/1.e3),'shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.tmin',str(tmin),'.pdf'])

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


text50_bothfree_pref 		= ''.join(['text_files/Pref/',name1_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedSP_pref  		= ''.join(['text_files/Pref/',name1_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text50_fixedOP_pref  		= ''.join(['text_files/Pref/',name1_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_bothfree_pref    	= ''.join(['text_files/Pref/',name3_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedSP_pref     	= ''.join(['text_files/Pref/',name3_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text250_fixedOP_pref		= ''.join(['text_files/Pref/',name3_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_bothfree_pref 		= ''.join(['text_files/Pref/',name4_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedSP_pref  		= ''.join(['text_files/Pref/',name4_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text500_fixedOP_pref 		= ''.join(['text_files/Pref/',name4_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_bothfree_pref		= ''.join(['text_files/Pref/',name5_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedSP_pref  		= ''.join(['text_files/Pref/',name5_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1000_fixedOP_pref  		= ''.join(['text_files/Pref/',name5_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_bothfree_pref    	= ''.join(['text_files/Pref/',name7_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedSP_pref     	= ''.join(['text_files/Pref/',name7_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text375_fixedOP_pref    	= ''.join(['text_files/Pref/',name7_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])


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
m50_bothfree_pref 	= np.loadtxt((text50_bothfree_pref))
m50_fixedSP_pref  	= np.loadtxt((text50_fixedSP_pref))
m50_fixedOP_pref  	= np.loadtxt((text50_fixedOP_pref))
m250_bothfree_pref 	= np.loadtxt((text250_bothfree_pref))
m250_fixedSP_pref 	= np.loadtxt((text250_fixedSP_pref))
m250_fixedOP_pref 	= np.loadtxt((text250_fixedOP_pref))           
m500_bothfree_pref 	= np.loadtxt((text500_bothfree_pref))
m500_fixedSP_pref  	= np.loadtxt((text500_fixedSP_pref))
m500_fixedOP_pref  	= np.loadtxt((text500_fixedOP_pref))    
m1000_bothfree_pref = np.loadtxt((text1000_bothfree_pref))
m1000_fixedSP_pref  = np.loadtxt((text1000_fixedSP_pref))
m1000_fixedOP_pref  = np.loadtxt((text1000_fixedOP_pref))    
m375_bothfree_pref 	= np.loadtxt((text375_bothfree_pref))
m375_fixedSP_pref 	= np.loadtxt((text375_fixedSP_pref))
m375_fixedOP_pref 	= np.loadtxt((text375_fixedOP_pref))

K_ind    = 11
ss_ind   = 17
sn_ind   = 6
anal_ind = 4
DP_ind   = 3
dip_ind  = 5
vc_ind   = 19
vs_ind   = 20

s_to_yr = 1./(365.25*24*3600)
cmyr_to_ms = 0.01/(365.25*24*3600)

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

fig=plt.figure()

# plot model full stress misfit vs. x-axis variable


# plot model DP misfit vs. x-axis variable
ax=fig.add_subplot(gs[0,0])

DP_array = np.empty((0, 2))

DP_array = append_DPcomponents(tmin,m50_fixedSP,m50_fixedOP_pref,DP_array)
DP_array = append_DPcomponents(tmin,m50_bothfree,m50_bothfree_pref,DP_array)
DP_array = append_DPcomponents(tmin,m50_fixedOP,m50_fixedOP_pref,DP_array)

DP_array = append_DPcomponents(tmin,m250_fixedSP,m250_fixedOP_pref,DP_array)
DP_array = append_DPcomponents(tmin,m250_bothfree,m250_bothfree_pref,DP_array)
DP_array = append_DPcomponents(tmin,m250_fixedOP,m250_fixedOP_pref,DP_array)

DP_array = append_DPcomponents(tmin,m375_fixedSP,m375_fixedOP_pref,DP_array)
DP_array = append_DPcomponents(tmin,m375_bothfree,m375_bothfree_pref,DP_array)
DP_array = append_DPcomponents(tmin,m375_fixedOP,m375_fixedOP_pref,DP_array)

DP_array = append_DPcomponents(tmin,m500_fixedSP,m500_fixedOP_pref,DP_array)
DP_array = append_DPcomponents(tmin,m500_bothfree,m500_bothfree_pref,DP_array)
DP_array = append_DPcomponents(tmin,m500_fixedOP,m500_fixedOP_pref,DP_array)  

DP_array = append_DPcomponents(tmin,m1000_fixedSP,m1000_fixedOP_pref,DP_array)
DP_array = append_DPcomponents(tmin,m1000_bothfree,m1000_bothfree_pref,DP_array)
DP_array = append_DPcomponents(tmin,m1000_fixedOP,m1000_fixedOP_pref ,DP_array)   



second_col = DP_array[:, 1]
second_col = second_col[(second_col >= -0.1) & (second_col <= 1.1)]
mean_val = np.mean(np.abs(second_col))
std_val = np.std(second_col)

bin_min = np.floor(second_col.min())
bin_max = np.ceil(second_col.max())
bins = np.arange(bin_min, bin_max + 0.1, 0.1)
plt.hist(second_col, bins=bins, edgecolor='black',zorder=2)

#annotate mean and sd
textstr = f'Mean = {mean_val:.2f}\nSD = {std_val:.2f}'
plt.annotate(textstr, xy=(0.95, 0.95), xycoords='axes fraction',
             fontsize=8, ha='left', va='top')

# axis stuff
#plt.ylim(0,  1); 
plt.xlim(-0.1,  1.1); 
plt.xlabel(r'$P_{subslab} / \Delta P$   [MPa]')
plt.ylabel(r'n')
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
plt.grid(True, which='both', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
plt.axvline(x=0.5, color='lightgray',linestyle='-',linewidth=1, zorder=0)
fixed_aspect_ratio(1)

plt.subplots_adjust(wspace=0.5)
plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=600)
plt.savefig(plot_name_pdf, format='pdf', bbox_inches='tight')
