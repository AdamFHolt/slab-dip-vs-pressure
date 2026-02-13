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
from functions import get_misfit_mean_and_stdev

analysis_depth_dz = float(sys.argv[1])     # m (depth for DP extraction and central point of shear stress derivative)
ds = float(sys.argv[2])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[3])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

analysis_depth1 = 230e3
analysis_depth2 = 330e3
analysis_depth3 = 430e3

plot_name_png = ''.join(['plots/DP-comparisons/compilations/averaged-misfit-vs-viscosity_all-mods.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/averaged-misfit-vs-viscosity_all-mods.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.pdf'])

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

# shallow
text1_shall_bothfree 	= ''.join(['text_files/',name1_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_shall_fixedSP  	= ''.join(['text_files/',name1_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_shall_fixedOP  	= ''.join(['text_files/',name1_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_shall_bothfree	= ''.join(['text_files/',name2_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_shall_fixedSP  	= ''.join(['text_files/',name2_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_shall_fixedOP  	= ''.join(['text_files/',name2_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_shall_bothfree    = ''.join(['text_files/',name3_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_shall_fixedSP     = ''.join(['text_files/',name3_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_shall_fixedOP	    = ''.join(['text_files/',name3_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_shall_bothfree 	= ''.join(['text_files/',name4_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_shall_fixedSP  	= ''.join(['text_files/',name4_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_shall_fixedOP  	= ''.join(['text_files/',name4_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_shall_bothfree	= ''.join(['text_files/',name5_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_shall_fixedSP  	= ''.join(['text_files/',name5_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_shall_fixedOP  	= ''.join(['text_files/',name5_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_shall_bothfree    = ''.join(['text_files/',name6_bothfree,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_shall_fixedSP     = ''.join(['text_files/',name6_fixedSP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_shall_fixedOP	    = ''.join(['text_files/',name6_fixedOP,'.z',str(analysis_depth1/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# intermediate
text1_int_bothfree 		= ''.join(['text_files/',name1_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_int_fixedSP  		= ''.join(['text_files/',name1_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_int_fixedOP  		= ''.join(['text_files/',name1_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_int_bothfree		= ''.join(['text_files/',name2_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_int_fixedSP  		= ''.join(['text_files/',name2_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_int_fixedOP  		= ''.join(['text_files/',name2_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_int_bothfree    	= ''.join(['text_files/',name3_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_int_fixedSP     	= ''.join(['text_files/',name3_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_int_fixedOP		= ''.join(['text_files/',name3_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_int_bothfree 		= ''.join(['text_files/',name4_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_int_fixedSP  		= ''.join(['text_files/',name4_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_int_fixedOP  		= ''.join(['text_files/',name4_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_int_bothfree		= ''.join(['text_files/',name5_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_int_fixedSP  		= ''.join(['text_files/',name5_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_int_fixedOP  		= ''.join(['text_files/',name5_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_int_bothfree    	= ''.join(['text_files/',name6_bothfree,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_int_fixedSP     	= ''.join(['text_files/',name6_fixedSP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_int_fixedOP	    = ''.join(['text_files/',name6_fixedOP,'.z',str(analysis_depth2/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
# deep
text1_deep_bothfree 	= ''.join(['text_files/',name1_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_deep_fixedSP  	= ''.join(['text_files/',name1_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text1_deep_fixedOP  	= ''.join(['text_files/',name1_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_deep_bothfree		= ''.join(['text_files/',name2_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_deep_fixedSP  	= ''.join(['text_files/',name2_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text2_deep_fixedOP  	= ''.join(['text_files/',name2_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_deep_bothfree    	= ''.join(['text_files/',name3_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_deep_fixedSP     	= ''.join(['text_files/',name3_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text3_deep_fixedOP	    = ''.join(['text_files/',name3_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_deep_bothfree 	= ''.join(['text_files/',name4_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_deep_fixedSP  	= ''.join(['text_files/',name4_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text4_deep_fixedOP  	= ''.join(['text_files/',name4_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_deep_bothfree		= ''.join(['text_files/',name5_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_deep_fixedSP  	= ''.join(['text_files/',name5_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text5_deep_fixedOP  	= ''.join(['text_files/',name5_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_deep_bothfree    	= ''.join(['text_files/',name6_bothfree,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_deep_fixedSP		= ''.join(['text_files/',name6_fixedSP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
text6_deep_fixedOP	    = ''.join(['text_files/',name6_fixedOP,'.z',str(analysis_depth3/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

# load in models
m1_shall_bothfree 	= np.loadtxt((text1_shall_bothfree)) 
m1_shall_fixedSP  	= np.loadtxt((text1_shall_fixedSP))
m1_shall_fixedOP  	= np.loadtxt((text1_shall_fixedOP))
m2_shall_bothfree 	= np.loadtxt((text2_shall_bothfree)) 
m2_shall_fixedSP  	= np.loadtxt((text2_shall_fixedSP))
m2_shall_fixedOP  	= np.loadtxt((text2_shall_fixedOP))
m3_shall_bothfree 	= np.loadtxt((text3_shall_bothfree)) 
m3_shall_fixedSP 	= np.loadtxt((text3_shall_fixedSP))
m3_shall_fixedOP 	= np.loadtxt((text3_shall_fixedOP))
m4_shall_bothfree 	= np.loadtxt((text4_shall_bothfree)) 
m4_shall_fixedSP  	= np.loadtxt((text4_shall_fixedSP))
m4_shall_fixedOP  	= np.loadtxt((text4_shall_fixedOP))
m5_shall_bothfree 	= np.loadtxt((text5_shall_bothfree)) 
m5_shall_fixedSP  	= np.loadtxt((text5_shall_fixedSP))
m5_shall_fixedOP  	= np.loadtxt((text5_shall_fixedOP))
m6_shall_bothfree 	= np.loadtxt((text6_shall_bothfree)) 
m6_shall_fixedSP 	= np.loadtxt((text6_shall_fixedSP))
m6_shall_fixedOP 	= np.loadtxt((text6_shall_fixedOP))
# --
m1_int_bothfree 	= np.loadtxt((text1_int_bothfree)) 
m1_int_fixedSP  	= np.loadtxt((text1_int_fixedSP))
m1_int_fixedOP  	= np.loadtxt((text1_int_fixedOP))
m2_int_bothfree 	= np.loadtxt((text2_int_bothfree)) 
m2_int_fixedSP  	= np.loadtxt((text2_int_fixedSP))
m2_int_fixedOP  	= np.loadtxt((text2_int_fixedOP))
m3_int_bothfree 	= np.loadtxt((text3_int_bothfree)) 
m3_int_fixedSP 		= np.loadtxt((text3_int_fixedSP))
m3_int_fixedOP 		= np.loadtxt((text3_int_fixedOP))
m4_int_bothfree 	= np.loadtxt((text4_int_bothfree)) 
m4_int_fixedSP  	= np.loadtxt((text4_int_fixedSP))
m4_int_fixedOP  	= np.loadtxt((text4_int_fixedOP))
m5_int_bothfree 	= np.loadtxt((text5_int_bothfree)) 
m5_int_fixedSP  	= np.loadtxt((text5_int_fixedSP))
m5_int_fixedOP  	= np.loadtxt((text5_int_fixedOP))
m6_int_bothfree 	= np.loadtxt((text6_int_bothfree)) 
m6_int_fixedSP 		= np.loadtxt((text6_int_fixedSP))
m6_int_fixedOP 		= np.loadtxt((text6_int_fixedOP))
# --
m1_deep_bothfree 	= np.loadtxt((text1_deep_bothfree)) 
m1_deep_fixedSP  	= np.loadtxt((text1_deep_fixedSP))
m1_deep_fixedOP  	= np.loadtxt((text1_deep_fixedOP))
m2_deep_bothfree 	= np.loadtxt((text2_deep_bothfree)) 
m2_deep_fixedSP  	= np.loadtxt((text2_deep_fixedSP))
m2_deep_fixedOP  	= np.loadtxt((text2_deep_fixedOP))
m3_deep_bothfree 	= np.loadtxt((text3_deep_bothfree)) 
m3_deep_fixedSP 	= np.loadtxt((text3_deep_fixedSP))
m3_deep_fixedOP 	= np.loadtxt((text3_deep_fixedOP))
m4_deep_bothfree 	= np.loadtxt((text4_deep_bothfree)) 
m4_deep_fixedSP  	= np.loadtxt((text4_deep_fixedSP))
m4_deep_fixedOP  	= np.loadtxt((text4_deep_fixedOP))
m5_deep_bothfree 	= np.loadtxt((text5_deep_bothfree)) 
m5_deep_fixedSP  	= np.loadtxt((text5_deep_fixedSP))
m5_deep_fixedOP  	= np.loadtxt((text5_deep_fixedOP))
m6_deep_bothfree 	= np.loadtxt((text6_deep_bothfree)) 
m6_deep_fixedSP 	= np.loadtxt((text6_deep_fixedSP))
m6_deep_fixedOP 	= np.loadtxt((text6_deep_fixedOP))
# --

# get misfit means and st. deviations
mean_woshear_m1_shall_bothfree, mean_wshear_m1_shall_bothfree, stdev_woshear_m1_shall_bothfree, stdev_wshear_m1_shall_bothfree 	= get_misfit_mean_and_stdev(m1_shall_bothfree)
mean_woshear_m1_shall_fixedSP,  mean_wshear_m1_shall_fixedSP,  stdev_woshear_m1_shall_fixedSP,  stdev_wshear_m1_shall_fixedSP 	= get_misfit_mean_and_stdev(m1_shall_fixedSP)
mean_woshear_m1_shall_fixedOP,  mean_wshear_m1_shall_fixedOP,  stdev_woshear_m1_shall_fixedOP,  stdev_wshear_m1_shall_fixedOP 	= get_misfit_mean_and_stdev(m1_shall_fixedOP)
mean_woshear_m2_shall_bothfree, mean_wshear_m2_shall_bothfree, stdev_woshear_m2_shall_bothfree, stdev_wshear_m2_shall_bothfree  = get_misfit_mean_and_stdev(m2_shall_bothfree)
mean_woshear_m2_shall_fixedSP,  mean_wshear_m2_shall_fixedSP,  stdev_woshear_m2_shall_fixedSP,  stdev_wshear_m2_shall_fixedSP 	= get_misfit_mean_and_stdev(m2_shall_fixedSP)
mean_woshear_m2_shall_fixedOP,  mean_wshear_m2_shall_fixedOP,  stdev_woshear_m2_shall_fixedOP,  stdev_wshear_m2_shall_fixedOP 	= get_misfit_mean_and_stdev(m2_shall_fixedOP)
mean_woshear_m3_shall_bothfree, mean_wshear_m3_shall_bothfree, stdev_woshear_m3_shall_bothfree, stdev_wshear_m3_shall_bothfree 	= get_misfit_mean_and_stdev(m3_shall_bothfree)
mean_woshear_m3_shall_fixedSP,  mean_wshear_m3_shall_fixedSP,  stdev_woshear_m3_shall_fixedSP,  stdev_wshear_m3_shall_fixedSP 	= get_misfit_mean_and_stdev(m3_shall_fixedSP)
mean_woshear_m3_shall_fixedOP,  mean_wshear_m3_shall_fixedOP,  stdev_woshear_m3_shall_fixedOP,  stdev_wshear_m3_shall_fixedOP 	= get_misfit_mean_and_stdev(m3_shall_fixedOP)
mean_woshear_m4_shall_bothfree, mean_wshear_m4_shall_bothfree, stdev_woshear_m4_shall_bothfree, stdev_wshear_m4_shall_bothfree 	= get_misfit_mean_and_stdev(m4_shall_bothfree)
mean_woshear_m4_shall_fixedSP,  mean_wshear_m4_shall_fixedSP,  stdev_woshear_m4_shall_fixedSP,  stdev_wshear_m4_shall_fixedSP 	= get_misfit_mean_and_stdev(m4_shall_fixedSP)
mean_woshear_m4_shall_fixedOP,  mean_wshear_m4_shall_fixedOP,  stdev_woshear_m4_shall_fixedOP,  stdev_wshear_m4_shall_fixedOP 	= get_misfit_mean_and_stdev(m4_shall_fixedOP)
mean_woshear_m5_shall_bothfree, mean_wshear_m5_shall_bothfree, stdev_woshear_m5_shall_bothfree, stdev_wshear_m5_shall_bothfree 	= get_misfit_mean_and_stdev(m5_shall_bothfree)
mean_woshear_m5_shall_fixedSP,  mean_wshear_m5_shall_fixedSP,  stdev_woshear_m5_shall_fixedSP,  stdev_wshear_m5_shall_fixedSP 	= get_misfit_mean_and_stdev(m5_shall_fixedSP)
mean_woshear_m5_shall_fixedOP,  mean_wshear_m5_shall_fixedOP,  stdev_woshear_m5_shall_fixedOP,  stdev_wshear_m5_shall_fixedOP 	= get_misfit_mean_and_stdev(m5_shall_fixedOP)
mean_woshear_m6_shall_bothfree, mean_wshear_m6_shall_bothfree, stdev_woshear_m6_shall_bothfree, stdev_wshear_m6_shall_bothfree 	= get_misfit_mean_and_stdev(m6_shall_bothfree)
mean_woshear_m6_shall_fixedSP,  mean_wshear_m6_shall_fixedSP,  stdev_woshear_m6_shall_fixedSP,  stdev_wshear_m6_shall_fixedSP 	= get_misfit_mean_and_stdev(m6_shall_fixedSP)
mean_woshear_m6_shall_fixedOP,  mean_wshear_m6_shall_fixedOP,  stdev_woshear_m6_shall_fixedOP,  stdev_wshear_m6_shall_fixedOP 	= get_misfit_mean_and_stdev(m6_shall_fixedOP)

mean_woshear_m1_int_bothfree, mean_wshear_m1_int_bothfree, stdev_woshear_m1_int_bothfree, stdev_wshear_m1_int_bothfree 	= get_misfit_mean_and_stdev(m1_int_bothfree)
mean_woshear_m1_int_fixedSP,  mean_wshear_m1_int_fixedSP,  stdev_woshear_m1_int_fixedSP,  stdev_wshear_m1_int_fixedSP 	= get_misfit_mean_and_stdev(m1_int_fixedSP)
mean_woshear_m1_int_fixedOP,  mean_wshear_m1_int_fixedOP,  stdev_woshear_m1_int_fixedOP,  stdev_wshear_m1_int_fixedOP 	= get_misfit_mean_and_stdev(m1_int_fixedOP)
mean_woshear_m2_int_bothfree, mean_wshear_m2_int_bothfree, stdev_woshear_m2_int_bothfree, stdev_wshear_m2_int_bothfree  = get_misfit_mean_and_stdev(m2_int_bothfree)
mean_woshear_m2_int_fixedSP,  mean_wshear_m2_int_fixedSP,  stdev_woshear_m2_int_fixedSP,  stdev_wshear_m2_int_fixedSP 	= get_misfit_mean_and_stdev(m2_int_fixedSP)
mean_woshear_m2_int_fixedOP,  mean_wshear_m2_int_fixedOP,  stdev_woshear_m2_int_fixedOP,  stdev_wshear_m2_int_fixedOP 	= get_misfit_mean_and_stdev(m2_int_fixedOP)
mean_woshear_m3_int_bothfree, mean_wshear_m3_int_bothfree, stdev_woshear_m3_int_bothfree, stdev_wshear_m3_int_bothfree 	= get_misfit_mean_and_stdev(m3_int_bothfree)
mean_woshear_m3_int_fixedSP,  mean_wshear_m3_int_fixedSP,  stdev_woshear_m3_int_fixedSP,  stdev_wshear_m3_int_fixedSP 	= get_misfit_mean_and_stdev(m3_int_fixedSP)
mean_woshear_m3_int_fixedOP,  mean_wshear_m3_int_fixedOP,  stdev_woshear_m3_int_fixedOP,  stdev_wshear_m3_int_fixedOP 	= get_misfit_mean_and_stdev(m3_int_fixedOP)
mean_woshear_m4_int_bothfree, mean_wshear_m4_int_bothfree, stdev_woshear_m4_int_bothfree, stdev_wshear_m4_int_bothfree 	= get_misfit_mean_and_stdev(m4_int_bothfree)
mean_woshear_m4_int_fixedSP,  mean_wshear_m4_int_fixedSP,  stdev_woshear_m4_int_fixedSP,  stdev_wshear_m4_int_fixedSP 	= get_misfit_mean_and_stdev(m4_int_fixedSP)
mean_woshear_m4_int_fixedOP,  mean_wshear_m4_int_fixedOP,  stdev_woshear_m4_int_fixedOP,  stdev_wshear_m4_int_fixedOP 	= get_misfit_mean_and_stdev(m4_int_fixedOP)
mean_woshear_m5_int_bothfree, mean_wshear_m5_int_bothfree, stdev_woshear_m5_int_bothfree, stdev_wshear_m5_int_bothfree 	= get_misfit_mean_and_stdev(m5_int_bothfree)
mean_woshear_m5_int_fixedSP,  mean_wshear_m5_int_fixedSP,  stdev_woshear_m5_int_fixedSP,  stdev_wshear_m5_int_fixedSP 	= get_misfit_mean_and_stdev(m5_int_fixedSP)
mean_woshear_m5_int_fixedOP,  mean_wshear_m5_int_fixedOP,  stdev_woshear_m5_int_fixedOP,  stdev_wshear_m5_int_fixedOP 	= get_misfit_mean_and_stdev(m5_int_fixedOP)
mean_woshear_m6_int_bothfree, mean_wshear_m6_int_bothfree, stdev_woshear_m6_int_bothfree, stdev_wshear_m6_int_bothfree 	= get_misfit_mean_and_stdev(m6_int_bothfree)
mean_woshear_m6_int_fixedSP,  mean_wshear_m6_int_fixedSP,  stdev_woshear_m6_int_fixedSP,  stdev_wshear_m6_int_fixedSP 	= get_misfit_mean_and_stdev(m6_int_fixedSP)
mean_woshear_m6_int_fixedOP,  mean_wshear_m6_int_fixedOP,  stdev_woshear_m6_int_fixedOP,  stdev_wshear_m6_int_fixedOP 	= get_misfit_mean_and_stdev(m6_int_fixedOP)

mean_woshear_m1_deep_bothfree, mean_wshear_m1_deep_bothfree, stdev_woshear_m1_deep_bothfree, stdev_wshear_m1_deep_bothfree 	= get_misfit_mean_and_stdev(m1_deep_bothfree)
mean_woshear_m1_deep_fixedSP,  mean_wshear_m1_deep_fixedSP,  stdev_woshear_m1_deep_fixedSP,  stdev_wshear_m1_deep_fixedSP 	= get_misfit_mean_and_stdev(m1_deep_fixedSP)
mean_woshear_m1_deep_fixedOP,  mean_wshear_m1_deep_fixedOP,  stdev_woshear_m1_deep_fixedOP,  stdev_wshear_m1_deep_fixedOP 	= get_misfit_mean_and_stdev(m1_deep_fixedOP)
mean_woshear_m2_deep_bothfree, mean_wshear_m2_deep_bothfree, stdev_woshear_m2_deep_bothfree, stdev_wshear_m2_deep_bothfree  = get_misfit_mean_and_stdev(m2_deep_bothfree)
mean_woshear_m2_deep_fixedSP,  mean_wshear_m2_deep_fixedSP,  stdev_woshear_m2_deep_fixedSP,  stdev_wshear_m2_deep_fixedSP 	= get_misfit_mean_and_stdev(m2_deep_fixedSP)
mean_woshear_m2_deep_fixedOP,  mean_wshear_m2_deep_fixedOP,  stdev_woshear_m2_deep_fixedOP,  stdev_wshear_m2_deep_fixedOP 	= get_misfit_mean_and_stdev(m2_deep_fixedOP)
mean_woshear_m3_deep_bothfree, mean_wshear_m3_deep_bothfree, stdev_woshear_m3_deep_bothfree, stdev_wshear_m3_deep_bothfree 	= get_misfit_mean_and_stdev(m3_deep_bothfree)
mean_woshear_m3_deep_fixedSP,  mean_wshear_m3_deep_fixedSP,  stdev_woshear_m3_deep_fixedSP,  stdev_wshear_m3_deep_fixedSP 	= get_misfit_mean_and_stdev(m3_deep_fixedSP)
mean_woshear_m3_deep_fixedOP,  mean_wshear_m3_deep_fixedOP,  stdev_woshear_m3_deep_fixedOP,  stdev_wshear_m3_deep_fixedOP 	= get_misfit_mean_and_stdev(m3_deep_fixedOP)
mean_woshear_m4_deep_bothfree, mean_wshear_m4_deep_bothfree, stdev_woshear_m4_deep_bothfree, stdev_wshear_m4_deep_bothfree 	= get_misfit_mean_and_stdev(m4_deep_bothfree)
mean_woshear_m4_deep_fixedSP,  mean_wshear_m4_deep_fixedSP,  stdev_woshear_m4_deep_fixedSP,  stdev_wshear_m4_deep_fixedSP 	= get_misfit_mean_and_stdev(m4_deep_fixedSP)
mean_woshear_m4_deep_fixedOP,  mean_wshear_m4_deep_fixedOP,  stdev_woshear_m4_deep_fixedOP,  stdev_wshear_m4_deep_fixedOP 	= get_misfit_mean_and_stdev(m4_deep_fixedOP)
mean_woshear_m5_deep_bothfree, mean_wshear_m5_deep_bothfree, stdev_woshear_m5_deep_bothfree, stdev_wshear_m5_deep_bothfree 	= get_misfit_mean_and_stdev(m5_deep_bothfree)
mean_woshear_m5_deep_fixedSP,  mean_wshear_m5_deep_fixedSP,  stdev_woshear_m5_deep_fixedSP,  stdev_wshear_m5_deep_fixedSP 	= get_misfit_mean_and_stdev(m5_deep_fixedSP)
mean_woshear_m5_deep_fixedOP,  mean_wshear_m5_deep_fixedOP,  stdev_woshear_m5_deep_fixedOP,  stdev_wshear_m5_deep_fixedOP 	= get_misfit_mean_and_stdev(m5_deep_fixedOP)
mean_woshear_m6_deep_bothfree, mean_wshear_m6_deep_bothfree, stdev_woshear_m6_deep_bothfree, stdev_wshear_m6_deep_bothfree 	= get_misfit_mean_and_stdev(m6_deep_bothfree)
mean_woshear_m6_deep_fixedSP,  mean_wshear_m6_deep_fixedSP,  stdev_woshear_m6_deep_fixedSP,  stdev_wshear_m6_deep_fixedSP 	= get_misfit_mean_and_stdev(m6_deep_fixedSP)
mean_woshear_m6_deep_fixedOP,  mean_wshear_m6_deep_fixedOP,  stdev_woshear_m6_deep_fixedOP,  stdev_wshear_m6_deep_fixedOP 	= get_misfit_mean_and_stdev(m6_deep_fixedOP)







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


ax=fig.add_subplot(gs[0,0])
plt.scatter(50,mean_woshear_m1_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(50,mean_woshear_m1_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(50,mean_woshear_m1_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(100,mean_woshear_m2_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(100,mean_woshear_m2_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(100,mean_woshear_m2_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(250,mean_woshear_m3_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(250,mean_woshear_m3_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(250,mean_woshear_m3_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(500,mean_woshear_m4_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(500,mean_woshear_m4_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(500,mean_woshear_m4_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(1000,mean_woshear_m5_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(1000,mean_woshear_m5_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(1000,mean_woshear_m5_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(2000,mean_woshear_m6_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(2000,mean_woshear_m6_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(2000,mean_woshear_m6_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab strength",size=6)
plt.ylabel("extracted-analytical [MPa]",size=6)
plt.xlim(-100,  2100); plt.ylim(-2,  26)
fixed_aspect_ratio(1)
plt.annotate('230 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


# misfit vs. curvature
ax=fig.add_subplot(gs[0,1])
plt.scatter(50,mean_woshear_m1_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(50,mean_woshear_m1_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(50,mean_woshear_m1_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(100,mean_woshear_m2_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(100,mean_woshear_m2_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(100,mean_woshear_m2_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(250,mean_woshear_m3_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(250,mean_woshear_m3_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(250,mean_woshear_m3_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(500,mean_woshear_m4_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(500,mean_woshear_m4_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(500,mean_woshear_m4_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(1000,mean_woshear_m5_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(1000,mean_woshear_m5_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(1000,mean_woshear_m5_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(2000,mean_woshear_m6_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(2000,mean_woshear_m6_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(2000,mean_woshear_m6_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab strength",size=6)
# plt.ylabel("extracted-analytical [MPa]",size=6)
plt.xlim(-100,  2100); plt.ylim(-2,  26)
fixed_aspect_ratio(1)
plt.annotate('330 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


# misfit vs. curvature
ax=fig.add_subplot(gs[0,2])
plt.scatter(50,mean_woshear_m1_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(50,mean_woshear_m1_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(50,mean_woshear_m1_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(100,mean_woshear_m2_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(100,mean_woshear_m2_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(100,mean_woshear_m2_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(250,mean_woshear_m3_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(250,mean_woshear_m3_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(250,mean_woshear_m3_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(500,mean_woshear_m4_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(500,mean_woshear_m4_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(500,mean_woshear_m4_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(1000,mean_woshear_m5_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(1000,mean_woshear_m5_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(1000,mean_woshear_m5_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(2000,mean_woshear_m6_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(2000,mean_woshear_m6_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(2000,mean_woshear_m6_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab strength",size=6)
# plt.ylabel("extracted-analytical [MPa]",size=6)
plt.xlim(-100,  2100); plt.ylim(-2,  26)
fixed_aspect_ratio(1)
plt.annotate('430 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')

#------

ax=fig.add_subplot(gs[1,0])
plt.scatter(50,mean_wshear_m1_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(50,mean_wshear_m1_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(50,mean_wshear_m1_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(100,mean_wshear_m2_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(100,mean_wshear_m2_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(100,mean_wshear_m2_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(250,mean_wshear_m3_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(250,mean_wshear_m3_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(250,mean_wshear_m3_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(500,mean_wshear_m4_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(500,mean_wshear_m4_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(500,mean_wshear_m4_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(1000,mean_wshear_m5_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(1000,mean_wshear_m5_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(1000,mean_wshear_m5_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(2000,mean_wshear_m6_shall_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(2000,mean_wshear_m6_shall_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(2000,mean_wshear_m6_shall_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab strength",size=6)
plt.ylabel("extracted+shear-analytical [MPa]",size=6)
plt.xlim(-100,  2100); plt.ylim(-2,  26)
fixed_aspect_ratio(1)
plt.annotate('230 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


# misfit vs. curvature
ax=fig.add_subplot(gs[1,1])
plt.scatter(50,mean_wshear_m1_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(50,mean_wshear_m1_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(50,mean_wshear_m1_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(100,mean_wshear_m2_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(100,mean_wshear_m2_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(100,mean_wshear_m2_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(250,mean_wshear_m3_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(250,mean_wshear_m3_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(250,mean_wshear_m3_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(500,mean_wshear_m4_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(500,mean_wshear_m4_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(500,mean_wshear_m4_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(1000,mean_wshear_m5_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(1000,mean_wshear_m5_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(1000,mean_wshear_m5_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(2000,mean_wshear_m6_int_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(2000,mean_wshear_m6_int_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(2000,mean_wshear_m6_int_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab strength",size=6)
# plt.ylabel("extracted+shear-analytical [MPa]",size=6)
plt.xlim(-100,  2100); plt.ylim(-2,  26)
fixed_aspect_ratio(1)
plt.annotate('330 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


# misfit vs. curvature
ax=fig.add_subplot(gs[1,2])
plt.scatter(50,mean_wshear_m1_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(50,mean_wshear_m1_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(50,mean_wshear_m1_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(100,mean_wshear_m2_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(100,mean_wshear_m2_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(100,mean_wshear_m2_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(250,mean_wshear_m3_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(250,mean_wshear_m3_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(250,mean_wshear_m3_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(500,mean_wshear_m4_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(500,mean_wshear_m4_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(500,mean_wshear_m4_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(1000,mean_wshear_m5_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(1000,mean_wshear_m5_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(1000,mean_wshear_m5_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 
plt.scatter(2000,mean_wshear_m6_deep_bothfree,s=12,color='black',edgecolor='black',zorder=3) 
plt.scatter(2000,mean_wshear_m6_deep_fixedSP,s=12,color='black',edgecolor='black',marker='v',zorder=3) 
plt.scatter(2000,mean_wshear_m6_deep_fixedOP,s=12,color='black',edgecolor='black',marker='^',zorder=3) 


plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab strength",size=6)
# plt.ylabel("extracted+shear-analytical [MPa]",size=6)
plt.xlim(-100,  2100); plt.ylim(-2,  26)
fixed_aspect_ratio(1)
plt.annotate('430 km depth', xy=(0.025,0.93), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=7,color='k')


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')


