#!/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
from numpy import trapz
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from scipy.ndimage import uniform_filter1d
from functions import create_grid, get_slab_midplane, get_dip_slab_midplane
from functions import extract_horiz_prof, get_slablocation_from_horiz_prof
from functions import get_nearslab_stresses, get_dip_at_certain_depth
from functions import get_stress_profile, convert_to_slabnorm_shearstress
from functions import convert_to_slabnorm_normstress, convert_to_slabnorm_normstress_TEST
from functions import get_curvature_slab_midplane, get_curvature_at_certain_depth
from functions import get_platevels_from_horiz_prof, get_slabvisc_from_horiz_prof
from functions import get_stress_at_certain_depth, get_farfieldP

model_name=str(sys.argv[1])            
max_time=int(sys.argv[2])   # largest number in csv_outputs/ filenames
analysis_depth = float(sys.argv[3])     # m (depth for DP extraction and central point of shear stress derivative)
analysis_depth_dz = float(sys.argv[4])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[5])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[6])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)

# model properties
xmax=5800.e3
ymax=1450.e3
drho = 50. # kg/m3

# ASPECT output 
csvs_loc =  'csv_outputs/'
models_loc =  'raw_outputs/'
stats_file = ''.join([models_loc,str(model_name),'/statistics'])
model_output_dt  = 50 # output dt as set in ASPECT .prm file (for getting the dimensional time)
num_header_lines = 16 # num header lines in stats_files (for getting the dimensional time)
saved_stresses_name = ''.join(['text_files/Pref/',model_name,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
print(saved_stresses_name)

# where to put the plots
plot_loc   = ''.join(['plots/evolution/',str(model_name)])
if not os.path.exists(plot_loc):
    os.mkdir(plot_loc)

# column numbers of the relevant properties in the .csv file. 
visc_col=26;    vx_col=0;           vy_col=1;
c_crust_col=23; c_ulith_col = 24;   c_llith_col = 25;
P_col = 29;     x_col=30;           y_col = 31;
sxx_col = 3;    syy_col = 7;        sxy_col = 4; 


first_time=8
saved_stresses = np.zeros(((max_time-first_time),2)) 
                                                     
ind = 0 

for time in range(first_time,max_time,1):

    # plot name
    csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])

    print("-------")
    print("t = %.0f" % (time))
    model_data  = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

    y_center = 300.
    Pref = get_farfieldP(y_center,ymax,dz,model_data,P_col,x_col,y_col)

    saved_stresses[ind,:] = time, Pref
    
    ind = ind + 1

np.savetxt(saved_stresses_name, saved_stresses)
