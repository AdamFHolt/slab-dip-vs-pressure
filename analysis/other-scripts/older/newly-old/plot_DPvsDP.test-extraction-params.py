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
from functions_plotting import plot_BvsFullForce_wKthresh, plot_BvsDP_wKthresh

curvature_thresh = float(sys.argv[1])   # 0.00075 # 1/km

tactual_min = 11 # first time step to use
tmin = tactual_min - 8


analysis_depth = 300e3
plot_name_png = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.z300km.test-extraction-params.png'])
plot_name_pdf = ''.join(['plots/DP-comparisons/compilations/DP-vs-DP.z300km.test-extraction-params.pdf'])
prof_dz = 1.e3

name_weak_bothfree   = "2D_compositional_subd_lower-res_new_50plates"
name_weak_fixedSP    = "2D_compositional_subd_lower-res_new_FixedSP_50plates"
name_weak_fixedOP    = "2D_compositional_subd_lower-res_new_FixedOP_50plates"
name_ref_bothfree    = "2D_compositional_subd_lower-res_new_250plates"
name_ref_fixedSP     = "2D_compositional_subd_lower-res_new_FixedSP_250plates"
name_ref_fixedOP     = "2D_compositional_subd_lower-res_new_FixedOP_250plates"
name_ref2_bothfree   = "2D_compositional_subd_lower-res_new2"         # 500
name_ref2_fixedSP    = "2D_compositional_subd_FixedSP_lower-res_new2" # 500 
name_ref2_fixedOP    = "2D_compositional_subd_FixedOP_lower-res_new" # 500
name_strong_bothfree = "2D_compositional_subd_lower-res_new_1000plates"
name_strong_fixedSP  = "2D_compositional_subd_lower-res_new_FixedSP_1000plates2"
name_strong_fixedOP  = "2D_compositional_subd_lower-res_new_FixedOP_1000plates"
name_new_bothfree    = "2D_compositional_subd_lower-res_new_375plates"
name_new_fixedSP     = "2D_compositional_subd_lower-res_new_FixedSP_375plates"
name_new_fixedOP     = "2D_compositional_subd_lower-res_new_FixedOP_375plates"


# ref
analysis_depth_dz=10.e3; ds=10.e3; 
text1_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text1_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
# 2
analysis_depth_dz2=7.5e3; ds=10.e3; 
text2_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text2_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz2/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
# 3
analysis_depth_dz3=15.e3; ds=10.e3; 
text3_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text3_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz3/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
# 4
analysis_depth_dz=10.e3; ds2=7.5e3; 
text4_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text4_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds2/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
# 5
analysis_depth_dz5=5.e3; ds=10e3; 
text5_ref_bothfree  = ''.join(['text_files/',name_ref_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_ref_fixedSP   = ''.join(['text_files/',name_ref_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_ref_fixedOP   = ''.join(['text_files/',name_ref_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_ref2_bothfree  = ''.join(['text_files/',name_ref2_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_ref2_fixedSP   = ''.join(['text_files/',name_ref2_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_ref2_fixedOP   = ''.join(['text_files/',name_ref2_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_weak_bothfree = ''.join(['text_files/',name_weak_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_weak_fixedSP      = ''.join(['text_files/',name_weak_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_weak_fixedOP      = ''.join(['text_files/',name_weak_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_strong_bothfree= ''.join(['text_files/',name_strong_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_strong_fixedSP = ''.join(['text_files/',name_strong_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_strong_fixedOP    = ''.join(['text_files/',name_strong_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_new_bothfree=''.join(['text_files/',name_new_bothfree,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_new_fixedSP= ''.join(['text_files/',name_new_fixedSP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])
text5_new_fixedOP= ''.join(['text_files/',name_new_fixedOP,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz5/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(prof_dz/1.e3),'km.txt'])


ref1_bothfree   = np.loadtxt((text1_ref_bothfree)) 
ref1_fixedSP    = np.loadtxt((text1_ref_fixedSP))
ref1_fixedOP    = np.loadtxt((text1_ref_fixedOP))
ref1b_bothfree   = np.loadtxt((text1_ref2_bothfree)) 
ref1b_fixedSP    = np.loadtxt((text1_ref2_fixedSP))
ref1b_fixedOP    = np.loadtxt((text1_ref2_fixedOP))
weak1_bothfree  = np.loadtxt((text1_weak_bothfree)) 
weak1_fixedSP   = np.loadtxt((text1_weak_fixedSP))
weak1_fixedOP   = np.loadtxt((text1_weak_fixedOP))
strong1_bothfree = np.loadtxt((text1_strong_bothfree)) 
strong1_fixedSP     = np.loadtxt((text1_strong_fixedSP))
strong1_fixedOP     = np.loadtxt((text1_strong_fixedOP))
new1_bothfree = np.loadtxt((text1_new_bothfree)) 
new1_fixedSP     = np.loadtxt((text1_new_fixedSP))
new1_fixedOP     = np.loadtxt((text1_new_fixedOP))

ref2_bothfree   = np.loadtxt((text2_ref_bothfree)) 
ref2_fixedSP    = np.loadtxt((text2_ref_fixedSP))
ref2_fixedOP    = np.loadtxt((text2_ref_fixedOP))
ref2b_bothfree   = np.loadtxt((text2_ref2_bothfree)) 
ref2b_fixedSP    = np.loadtxt((text2_ref2_fixedSP))
ref2b_fixedOP    = np.loadtxt((text2_ref2_fixedOP))
weak2_bothfree  = np.loadtxt((text2_weak_bothfree)) 
weak2_fixedSP   = np.loadtxt((text2_weak_fixedSP))
weak2_fixedOP   = np.loadtxt((text2_weak_fixedOP))
strong2_bothfree = np.loadtxt((text2_strong_bothfree)) 
strong2_fixedSP     = np.loadtxt((text2_strong_fixedSP))
strong2_fixedOP     = np.loadtxt((text2_strong_fixedOP))
new2_bothfree = np.loadtxt((text2_new_bothfree)) 
new2_fixedSP     = np.loadtxt((text2_new_fixedSP))
new2_fixedOP     = np.loadtxt((text2_new_fixedOP))

ref3_bothfree   = np.loadtxt((text3_ref_bothfree)) 
ref3_fixedSP    = np.loadtxt((text3_ref_fixedSP))
ref3_fixedOP    = np.loadtxt((text3_ref_fixedOP))
ref3b_bothfree   = np.loadtxt((text3_ref2_bothfree)) 
ref3b_fixedSP    = np.loadtxt((text3_ref2_fixedSP))
ref3b_fixedOP    = np.loadtxt((text3_ref2_fixedOP))
weak3_bothfree  = np.loadtxt((text3_weak_bothfree)) 
weak3_fixedSP   = np.loadtxt((text3_weak_fixedSP))
weak3_fixedOP   = np.loadtxt((text3_weak_fixedOP))
strong3_bothfree = np.loadtxt((text3_strong_bothfree)) 
strong3_fixedSP     = np.loadtxt((text3_strong_fixedSP))
strong3_fixedOP     = np.loadtxt((text3_strong_fixedOP))
new3_bothfree = np.loadtxt((text3_new_bothfree)) 
new3_fixedSP     = np.loadtxt((text3_new_fixedSP))
new3_fixedOP     = np.loadtxt((text3_new_fixedOP))

ref4_bothfree   = np.loadtxt((text4_ref_bothfree)) 
ref4_fixedSP    = np.loadtxt((text4_ref_fixedSP))
ref4_fixedOP    = np.loadtxt((text4_ref_fixedOP))
ref4b_bothfree   = np.loadtxt((text4_ref2_bothfree)) 
ref4b_fixedSP    = np.loadtxt((text4_ref2_fixedSP))
ref4b_fixedOP    = np.loadtxt((text4_ref2_fixedOP))
weak4_bothfree  = np.loadtxt((text4_weak_bothfree)) 
weak4_fixedSP   = np.loadtxt((text4_weak_fixedSP))
weak4_fixedOP   = np.loadtxt((text4_weak_fixedOP))
strong4_bothfree = np.loadtxt((text4_strong_bothfree)) 
strong4_fixedSP     = np.loadtxt((text4_strong_fixedSP))
strong4_fixedOP     = np.loadtxt((text4_strong_fixedOP))
new4_bothfree = np.loadtxt((text4_new_bothfree)) 
new4_fixedSP     = np.loadtxt((text4_new_fixedSP))
new4_fixedOP     = np.loadtxt((text4_new_fixedOP))

ref5_bothfree   = np.loadtxt((text5_ref_bothfree)) 
ref5_fixedSP    = np.loadtxt((text5_ref_fixedSP))
ref5_fixedOP    = np.loadtxt((text5_ref_fixedOP))
ref5b_bothfree   = np.loadtxt((text5_ref2_bothfree)) 
ref5b_fixedSP    = np.loadtxt((text5_ref2_fixedSP))
ref5b_fixedOP    = np.loadtxt((text5_ref2_fixedOP))
weak5_bothfree  = np.loadtxt((text5_weak_bothfree)) 
weak5_fixedSP   = np.loadtxt((text5_weak_fixedSP))
weak5_fixedOP   = np.loadtxt((text5_weak_fixedOP))
strong5_bothfree = np.loadtxt((text5_strong_bothfree)) 
strong5_fixedSP     = np.loadtxt((text5_strong_fixedSP))
strong5_fixedOP     = np.loadtxt((text5_strong_fixedOP))
new5_bothfree = np.loadtxt((text5_new_bothfree)) 
new5_fixedSP     = np.loadtxt((text5_new_fixedSP))
new5_fixedOP     = np.loadtxt((text5_new_fixedOP))

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
# 50
plot_BvsFullForce_wKthresh(tmin,weak1_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak1_fixedSP,curvature_thresh, 'tan','black','v')
plot_BvsFullForce_wKthresh(tmin,weak1_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref1_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref1_fixedSP,curvature_thresh, 'peru','black','v')
plot_BvsFullForce_wKthresh(tmin,ref1_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new1_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new1_fixedSP,curvature_thresh, 'firebrick','black','v')
plot_BvsFullForce_wKthresh(tmin,new1_fixedOP,curvature_thresh, 'firebrick','black','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref1b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref1b_fixedSP,curvature_thresh, 'maroon','black','v')
plot_BvsFullForce_wKthresh(tmin,ref1b_fixedOP,curvature_thresh, 'maroon','black','^')
# 1000
plot_BvsFullForce_wKthresh(tmin,strong1_bothfree,curvature_thresh,'black','black','o')
plot_BvsFullForce_wKthresh(tmin,strong1_fixedSP,curvature_thresh, 'black','black','v')
plot_BvsFullForce_wKthresh(tmin,strong1_fixedOP,curvature_thresh, 'black','black','^')

# text
ax.annotate('dz=10 km, DP ds=10 km',size=6,
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')

# axis stuff
plt.xlim(0,  30); plt.ylim(0,  30)
plt.plot([0, 30], [0, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.ylabel(r'$\mathregular{\Delta P +}$ slab stress  [MPa]',size=6.5)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [0,10,20,30] )
fixed_aspect_ratio(1)


ax=fig.add_subplot(gs[0,1])
# 50
plot_BvsFullForce_wKthresh(tmin,weak2_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak2_fixedSP,curvature_thresh, 'tan','black','v')
plot_BvsFullForce_wKthresh(tmin,weak2_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref2_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref2_fixedSP,curvature_thresh, 'peru','black','v')
plot_BvsFullForce_wKthresh(tmin,ref2_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new2_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new2_fixedSP,curvature_thresh, 'firebrick','black','v')
plot_BvsFullForce_wKthresh(tmin,new2_fixedOP,curvature_thresh, 'firebrick','black','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref2b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref2b_fixedSP,curvature_thresh, 'maroon','black','v')
plot_BvsFullForce_wKthresh(tmin,ref2b_fixedOP,curvature_thresh, 'maroon','black','^')
# 1000
plot_BvsFullForce_wKthresh(tmin,strong2_bothfree,curvature_thresh,'black','black','o')
plot_BvsFullForce_wKthresh(tmin,strong2_fixedSP,curvature_thresh, 'black','black','v')
plot_BvsFullForce_wKthresh(tmin,strong2_fixedOP,curvature_thresh, 'black','black','^')

# text
ax.annotate('dz=7.5 km, DP ds=10 km',size=6,
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')

# axis stuff
plt.xlim(0,  30); plt.ylim(0,  30)
plt.plot([0, 30], [0, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [0,10,20,30] )
fixed_aspect_ratio(1)


ax=fig.add_subplot(gs[0,2])
# 50
plot_BvsFullForce_wKthresh(tmin,weak3_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak3_fixedSP,curvature_thresh, 'tan','black','v')
plot_BvsFullForce_wKthresh(tmin,weak3_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref3_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref3_fixedSP,curvature_thresh, 'peru','black','v')
plot_BvsFullForce_wKthresh(tmin,ref3_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new3_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new3_fixedSP,curvature_thresh, 'firebrick','black','v')
plot_BvsFullForce_wKthresh(tmin,new3_fixedOP,curvature_thresh, 'firebrick','black','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref3b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref3b_fixedSP,curvature_thresh, 'maroon','black','v')
plot_BvsFullForce_wKthresh(tmin,ref3b_fixedOP,curvature_thresh, 'maroon','black','^')
# 1000
plot_BvsFullForce_wKthresh(tmin,strong3_bothfree,curvature_thresh,'black','black','o')
plot_BvsFullForce_wKthresh(tmin,strong3_fixedSP,curvature_thresh, 'black','black','v')
plot_BvsFullForce_wKthresh(tmin,strong3_fixedOP,curvature_thresh, 'black','black','^')

# text
ax.annotate('dz=15 km, DP ds=10 km',size=6,
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')


# axis stuff
plt.xlim(0,  30); plt.ylim(0,  30)
plt.plot([0, 30], [0, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [0,10,20,30] )
fixed_aspect_ratio(1)

ax=fig.add_subplot(gs[1,0])
# 50
plot_BvsFullForce_wKthresh(tmin,weak4_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak4_fixedSP,curvature_thresh, 'tan','black','v')
plot_BvsFullForce_wKthresh(tmin,weak4_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref4_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref4_fixedSP,curvature_thresh, 'peru','black','v')
plot_BvsFullForce_wKthresh(tmin,ref4_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new4_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new4_fixedSP,curvature_thresh, 'firebrick','black','v')
plot_BvsFullForce_wKthresh(tmin,new4_fixedOP,curvature_thresh, 'firebrick','black','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref4b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref4b_fixedSP,curvature_thresh, 'maroon','black','v')
plot_BvsFullForce_wKthresh(tmin,ref4b_fixedOP,curvature_thresh, 'maroon','black','^')
# 1000
plot_BvsFullForce_wKthresh(tmin,strong4_bothfree,curvature_thresh,'black','black','o')
plot_BvsFullForce_wKthresh(tmin,strong4_fixedSP,curvature_thresh, 'black','black','v')
plot_BvsFullForce_wKthresh(tmin,strong4_fixedOP,curvature_thresh, 'black','black','^')

# text
ax.annotate('dz=10 km, DP ds=7.5 km',size=6,
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')

# axis stuff
plt.xlim(0,  30); plt.ylim(0,  30)
plt.plot([0, 30], [0, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [0,10,20,30] )
fixed_aspect_ratio(1)


ax=fig.add_subplot(gs[1,1])
# 50
plot_BvsFullForce_wKthresh(tmin,weak5_bothfree,curvature_thresh,'tan','black','o')
plot_BvsFullForce_wKthresh(tmin,weak5_fixedSP,curvature_thresh, 'tan','black','v')
plot_BvsFullForce_wKthresh(tmin,weak5_fixedOP,curvature_thresh, 'tan','black','^')
# 250
plot_BvsFullForce_wKthresh(tmin,ref5_bothfree,curvature_thresh,'peru','black','o')
plot_BvsFullForce_wKthresh(tmin,ref5_fixedSP,curvature_thresh, 'peru','black','v')
plot_BvsFullForce_wKthresh(tmin,ref5_fixedOP,curvature_thresh, 'peru','black','^')
# 375
plot_BvsFullForce_wKthresh(tmin,new5_bothfree,curvature_thresh,'firebrick','black','o')
plot_BvsFullForce_wKthresh(tmin,new5_fixedSP,curvature_thresh, 'firebrick','black','v')
plot_BvsFullForce_wKthresh(tmin,new5_fixedOP,curvature_thresh, 'firebrick','black','^')
# 500
plot_BvsFullForce_wKthresh(tmin,ref5b_bothfree,curvature_thresh,'maroon','black','o')
plot_BvsFullForce_wKthresh(tmin,ref5b_fixedSP,curvature_thresh, 'maroon','black','v')
plot_BvsFullForce_wKthresh(tmin,ref5b_fixedOP,curvature_thresh, 'maroon','black','^')
# 1000
plot_BvsFullForce_wKthresh(tmin,strong5_bothfree,curvature_thresh,'black','black','o')
plot_BvsFullForce_wKthresh(tmin,strong5_fixedSP,curvature_thresh, 'black','black','v')
plot_BvsFullForce_wKthresh(tmin,strong5_fixedOP,curvature_thresh, 'black','black','^')

# text
ax.annotate('dz=5 km, DP ds=10 km',size=6,
            xy=(1, 0), xycoords='axes fraction',
            xytext=(-20, 20), textcoords='offset pixels',
            horizontalalignment='right',
            verticalalignment='bottom')

# axis stuff
plt.xlim(0,  30); plt.ylim(0,  30)
plt.plot([0, 30], [0, 30], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel('slab buoyancy [MPa]',size=6.5)
ax.set_xticks( [0,10,20,30] )
fixed_aspect_ratio(1)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)
# plt.savefig(plot_name_pdf, bbox_inches='tight', format='pdf')



