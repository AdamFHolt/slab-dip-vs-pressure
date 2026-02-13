#!/bin/python 1

import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
from numpy import trapz
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev
from functions import create_grid, get_slab_midplane, get_dip_slab_midplane
from functions import extract_horiz_prof, get_slablocation_from_horiz_prof
from functions import get_nearslab_stresses, get_dip_at_certain_depth
from functions import get_stress_profile, convert_to_slabnorm_shearstress, convert_to_slabnorm_normstress

model_name=str(sys.argv[1])
analysis_depth = float(sys.argv[2])     #analysis_depth = 230.e3
analysis_depth_dz = 10.e3
ds = 7.5e3
dz = 1.e3

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
saved_stresses_name = ''.join(['text_files/lo-res/',model_name,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

# where to put the plots
plot_loc   = ''.join(['plots/tests/',str(model_name)])
if not os.path.exists(plot_loc):
	os.mkdir(plot_loc)

# column numbers of the relevant properties in the .csv file. 
visc_col=26;    vx_col=0;           vy_col=1;
c_ulith_col = 24;   c_llith_col = 25;
P_col = 29;     x_col=30;           y_col = 31;
sxx_col = 3;    syy_col = 7;        sxy_col = 4; 

# create low res grid for plotting
xmin_plot=0; ymin_plot=500.e3; grid_res=2.5e3
X_low, Y_low = create_grid(xmin_plot,xmax,ymin_plot,ymax,grid_res)

# higher res grid for detailed calculations
xmin_plot2=2500.e3; xmax_plot2=3750.e3
ymin_plot2=ymax-600.e3; grid_res2=0.25e3
X_low2, Y_low2 = create_grid(xmin_plot2,xmax_plot2,ymin_plot2,ymax,grid_res2)

ind = 0 

# for time in range(first_time,max_time+1,1):
times = np.array([14,16,18,20,22]) 

saved_stresses = np.zeros((len(times),12)) # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C,slabnorm thick, snorm_left, snorm_right

for time in times:

	# plot name
	csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])
	plotname=''.join([plot_loc,'/z',str(analysis_depth/1.e3),'km.',str(time),'.png'])

	# get dimensional time
	stats_line_num = num_header_lines + (time * model_output_dt) 
	f=open(stats_file); line=f.readlines()[stats_line_num]
	time_dim=float(line.split()[1])/1.e6 # Myr
	print("-------")
	print("%.0f: t = %.1f Myr" % (time,time_dim))
	model_data  = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

	# print "interpolating model outputs to regular grid..."
	visc   = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,visc_col], (X_low, Y_low), method='linear')
	llith  = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,c_llith_col], (X_low2, Y_low2), method='linear')

	# get lithosphere contour, and then trim it to get slab mid-plane
	comp_contour_val = 0.5
	llith_cont = plt.contour(X_low2/1.e3, (ymax-Y_low2)/1.e3, llith, levels=[comp_contour_val])
	llith_points_tmp = llith_cont.collections[0].get_paths()[0].vertices
	cutoff_shall = 110.; cutoff_deep  = 575.;
	llith_points = get_slab_midplane(llith_points_tmp,cutoff_shall,cutoff_deep)

	# get horizontal profiles and slab locations 
	midmant_prof = extract_horiz_prof(model_data,analysis_depth,ymax,dz,x_col,y_col)
	x_left, y_left, x_center, y_center, x_right, y_right = \
		get_slablocation_from_horiz_prof(midmant_prof,x_col,y_col,c_ulith_col,c_llith_col,ymax)

	analysis_depth_shall = analysis_depth - analysis_depth_dz
	midmant_prof_shall   = extract_horiz_prof(model_data,analysis_depth_shall,ymax,dz,x_col,y_col)
	x_left_shall, y_left_shall, x_center_shall, y_center_shall, x_right_shall, y_right_shall = \
		get_slablocation_from_horiz_prof(midmant_prof_shall,x_col,y_col,c_ulith_col,c_llith_col,ymax)
		
	analysis_depth_deep  = analysis_depth + analysis_depth_dz
	midmant_prof_deep    = extract_horiz_prof(model_data,analysis_depth_deep,ymax,dz,x_col,y_col)
	x_left_deep,  y_left_deep,  x_center_deep,  y_center_deep,  x_right_deep,  y_right_deep = \
		get_slablocation_from_horiz_prof(midmant_prof_deep,x_col,y_col,c_ulith_col,c_llith_col,ymax)

	# get dip of mid-plane (and smooth)
	dips = get_dip_slab_midplane(llith_points)
	# dips_filtered = savgol_filter(dips[:,0],151,3)
	dips = savgol_filter(dips[:,0],901,3)

	# get dips at the relevant depths
	dip_midmant, xloc_dip             = get_dip_at_certain_depth(dips,llith_points,analysis_depth)
	dip_midmant_shall, xloc_dip_shall = get_dip_at_certain_depth(dips,llith_points,analysis_depth_shall)
	dip_midmant_deep, xloc_dip_deep   = get_dip_at_certain_depth(dips,llith_points,analysis_depth_deep)
	print("dips: shallow = %.2f, center = %.2f, deep = %.2f deg" % (dip_midmant_shall,dip_midmant,dip_midmant_deep))


	# slab normal thicknesses at the relevant depths
	slabnorm_thick       = (x_right-x_left)             * np.sin(np.deg2rad(dip_midmant))       # m
	slabnorm_thick_ll    = (x_center-x_left)             * np.sin(np.deg2rad(dip_midmant))       # m
	slabnorm_thick_ul    = (x_right-x_center)             * np.sin(np.deg2rad(dip_midmant))       # m
	# print(slabnorm_thick,slabnorm_thick_ll,slabnorm_thick_ul)

	slabnorm_thick_shall    = (x_right_shall-x_left_shall) * np.sin(np.deg2rad(dip_midmant_shall)) # m
	slabnorm_thick_shall_ll = (x_center_shall-x_left_shall) * np.sin(np.deg2rad(dip_midmant_shall)) # m
	slabnorm_thick_shall_ul = (x_right_shall-x_center_shall) * np.sin(np.deg2rad(dip_midmant_shall)) # m
	# print(slabnorm_thick_shall,slabnorm_thick_shall_ll,slabnorm_thick_shall_ul)

	slabnorm_thick_deep     = (x_right_deep -x_left_deep)  * np.sin(np.deg2rad(dip_midmant_deep))  # m
	slabnorm_thick_deep_ll  = (x_center_deep -x_left_deep)  * np.sin(np.deg2rad(dip_midmant_deep))  # m
	slabnorm_thick_deep_ul  = (x_right_deep -x_center_deep)  * np.sin(np.deg2rad(dip_midmant_deep))  # m
	# print(slabnorm_thick_deep,slabnorm_thick_deep_ll,slabnorm_thick_deep_ul)

	# compute theoretical DP
	DP_anal = drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(dip_midmant)) # Pa

	# pull out model DP
	Pleft, Pright, Sxxleft, Sxxright, Sxyleft, Sxyright, Pleft_x, Pright_x, Pleft_y, Pright_y, slab_x_left, slab_y_left, slab_x_right, slab_y_right = \
		get_nearslab_stresses(y_center,slabnorm_thick_ul,slabnorm_thick_ll,dip_midmant,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz) # Pa
	DP_mod = Pleft - Pright # Pa
	DP_missing = DP_anal - DP_mod
	print("DPs: analytical = %.1f, model %.1f, missing = %.1f MPa" % (DP_anal/1.e6,DP_mod/1.e6,DP_missing/1.e6))
	Snorm_left  = Sxyleft*np.sin(np.deg2rad(2.*dip_midmant))  + Sxxleft*((np.sin(np.deg2rad(dip_midmant)))**2 - (np.cos(np.deg2rad(dip_midmant)))**2)
	Snorm_right = Sxyright*np.sin(np.deg2rad(2.*dip_midmant)) + Sxxright*((np.sin(np.deg2rad(dip_midmant)))**2 - (np.cos(np.deg2rad(dip_midmant)))**2)
	Snorm_contrib = Snorm_left-Snorm_right # positive = upwards push (helping to support buoyancy) 
	print("Dev stresses at the slab surfaces: %.2f MPa upwards push" % (Snorm_contrib/1.e6))

	# other depths for shear stress
	Pleft_shall, Pright_shall, Sxxleft_shall, Sxxright_shall, Sxyleft_shall, Sxyright_shall, Pleft_x_shall, Pright_x_shall, Pleft_y_shall, Pright_y_shall, slab_x_left_shall, slab_y_left_shall, slab_x_right_shall, slab_y_right_shall = \
		get_nearslab_stresses(y_center_shall,slabnorm_thick_ul,slabnorm_thick_ll,dip_midmant_shall,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz) # Pa
	Pleft_deep, Pright_deep, Sxxleft_deep, Sxxright_deep, Sxyleft_deep, Sxyright_deep, Pleft_x_deep, Pright_x_deep, Pleft_y_deep, Pright_y_deep, slab_x_left_deep, slab_y_left_deep, slab_x_right_deep, slab_y_right_deep = \
		get_nearslab_stresses(y_center_deep, slabnorm_thick_ul,slabnorm_thick_ll, dip_midmant_deep, model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz) # Pa  
	DP_mod_shall = Pleft_shall - Pright_shall; DP_mod_deep  = Pleft_deep  - Pright_deep

	# get slab stresses along slab-perp profiles
	sxx_slab, sxy_slab, P_slab, profile = \
		get_stress_profile(model_data,slab_x_left,slab_y_left,slab_x_right,slab_y_right,ymax,x_col,y_col,sxx_col,sxy_col,P_col)
	sxx_slab_shall, sxy_slab_shall, P_slab_shall, profile_shall = \
		get_stress_profile(model_data,slab_x_left_shall,slab_y_left_shall,slab_x_right_shall,slab_y_right_shall,ymax,x_col,y_col,sxx_col,sxy_col,P_col)
	sxx_slab_deep, sxy_slab_deep, P_slab_deep, profile_deep = \
		get_stress_profile(model_data,slab_x_left_deep,slab_y_left_deep,slab_x_right_deep,slab_y_right_deep,ymax,x_col,y_col,sxx_col,sxy_col,P_col)

	# rotate to get shear stresses
	slab_shear_stress       = convert_to_slabnorm_shearstress(sxx_slab,      sxy_slab,      dip_midmant)
	slab_shear_stress_shall = convert_to_slabnorm_shearstress(sxx_slab_shall,sxy_slab_shall,dip_midmant_shall)
	slab_shear_stress_deep  = convert_to_slabnorm_shearstress(sxx_slab_deep, sxy_slab_deep, dip_midmant_deep)

	# rotate to get slab-normal normal stress
	slab_norm_devstress,       slab_norm_fullstress    		= convert_to_slabnorm_normstress(sxx_slab, sxy_slab, P_slab,  dip_midmant)
	slab_norm_devstress_shall, slab_norm_fullstress_shall 	= convert_to_slabnorm_normstress(sxx_slab_shall, sxy_slab_shall, P_slab_shall,  dip_midmant_shall)
	slab_norm_devstress_deep,  slab_norm_fullstress_deep 	= convert_to_slabnorm_normstress(sxx_slab_deep,  sxy_slab_deep,  P_slab_deep,   dip_midmant_deep)

	# cut out the boundaries of the slabs (which have stress discontinuities)
	profcut_min = 10 # 7
	profcut_max = 10 # 5
	inds                        = np.where((profile[:,2]>profcut_min) & (profile[:,2]<profile[len(profile)-1,2]-profcut_max))[0]
	slab_shear_stress_cut       = slab_shear_stress[inds]
	profile_cut                 = profile[inds,:]
	# inds_shall                  = np.where((profile_shall[:,2]>profcut_min) & (profile_shall[:,2]<profile_shall[len(profile_shall)-1,2]-profcut_max))[0]
	slab_shear_stress_shall_cut = slab_shear_stress_shall[inds]
	profile_shall_cut           = profile_shall[inds,:]
	# inds_deep                   = np.where((profile_deep[:,2]>profcut_min) & (profile_deep[:,2]<profile_deep[len(profile_deep)-1,2]-profcut_max))[0]
	slab_shear_stress_deep_cut  = slab_shear_stress_deep[inds]
	profile_deep_cut            = profile_deep[inds,:]

	# integrate stresses to get slab-perp shear force
	slab_force          = trapz(slab_shear_stress_cut[:,0],         profile_cut[:,2]*1.e3)
	slab_force_shall    = trapz(slab_shear_stress_shall_cut[:,0],   profile_shall_cut[:,2]*1.e3)
	slab_force_deep     = trapz(slab_shear_stress_deep_cut[:,0],    profile_deep_cut[:,2]*1.e3)

	# convert to stress 
	prof_length         = profile_cut[len(profile_cut)-1,2] - profile_cut[0,2]
	slab_stress         = slab_force / (prof_length*1.e3)
	prof_length_shall   = profile_shall_cut[len(profile_shall_cut)-1,2] - profile_shall_cut[0,2]
	slab_stress_shall   = slab_force_shall / (prof_length_shall*1.e3)
	prof_length_deep    = profile_deep_cut[len(profile_deep_cut)-1,2] - profile_deep_cut[0,2]
	slab_stress_deep    = slab_force_deep  / (prof_length_deep*1.e3)
	# print("shear stresses: shallow = %.1f; center = %.1f; deep = %.1f MPa" % (slab_stress_shall/1.e6,slab_stress/1.e6,slab_stress_deep/1.e6))
   
	ds_shall = np.sqrt((x_center - x_center_shall)**2 + (y_center*1e3 - y_center_shall*1e3)**2)
	ds_deep  = np.sqrt((x_center - x_center_deep)**2  + (y_center*1e3 - y_center_deep*1e3)**2)

	slab_stress_gradient = (slab_stress_deep - slab_stress_shall)/(ds_shall + ds_deep)
	slab_stress_term = slab_stress_gradient * 0.5*(slabnorm_thick_shall + slabnorm_thick_deep)
	slab_stress_gradient_b = (slab_stress_deep - slab_stress)/(ds_deep)
	slab_stress_term_b = slab_stress_gradient_b * 0.5*(slabnorm_thick + slabnorm_thick_deep)
	slab_stress_gradient_c = (slab_stress - slab_stress_shall)/(ds_shall)
	slab_stress_term_c = slab_stress_gradient_c * 0.5*(slabnorm_thick + slabnorm_thick_shall)

	# # alternative shear stress calculation
	# slab_stress_fullgradient = (slab_shear_stress_deep_cut-slab_shear_stress_shall_cut)/(ds_shall+ds_deep)
	# slab_stress_intgradient  = trapz(slab_stress_fullgradient[:,0], profile_cut[:,2]*1.e3)
	# avg_prof_length = 0.5*(prof_length_shall+prof_length_deep)*1.e3
	# avg_slab_thick  = 0.5*(slabnorm_thick_shall + slabnorm_thick_deep)
	# slab_stress_intgradient = slab_stress_intgradient * (avg_slab_thick/avg_prof_length)

	print("shear stress term = %.2f MPa (others: %.2f, %.2f)" % (slab_stress_term/1.e6,slab_stress_term_b/1.e6,slab_stress_term_c/1.e6))
	# print("shear stress term 2 = %.2f MPa" % (slab_stress_intgradient/1.e6))

	# normal stress gradient calculation:
	slab_norm_fullstress_cut 		= slab_norm_fullstress[inds]
	slab_norm_fullstress_shall_cut 	= slab_norm_fullstress_shall[inds]
	slab_norm_fullstress_deep_cut  	= slab_norm_fullstress_deep[inds]
	slab_norm_fullstress_gradient 		= np.zeros(len(slab_norm_fullstress_cut))
	slab_norm_fullstress_shall_gradient = np.zeros(len(slab_norm_fullstress_shall_cut))
	slab_norm_fullstress_deep_gradient 	= np.zeros(len(slab_norm_fullstress_deep_cut))
	dr = 1e9; dr_shall = 1e9; dr_deep = 1e9 
	gradient_min = 1e9; gradient_max = -1e9;
	gradient_deep_min = 1e9;  gradient_deep_max = -1e9;
	gradient_shall_min = 1e9; gradient_shall_max = -1e9;
	for j in range(len(slab_norm_fullstress_cut)):
		if j == 0:
			dstress = (slab_norm_fullstress_cut[j+1] - slab_norm_fullstress_cut[j])/1.e6
			dprof   = np.sqrt((profile_cut[j+1,0]-profile_cut[j,0])**2 + (profile_cut[j+1,1]-profile_cut[j,1])**2)
			slab_norm_fullstress_gradient[j] = dstress/dprof  # MPa / km
			# shall
			dstress_shall = (slab_norm_fullstress_shall_cut[j+1] - slab_norm_fullstress_shall_cut[j])/1.e6
			dprof_shall   = np.sqrt((profile_shall_cut[j+1,0]-profile_shall_cut[j,0])**2 + (profile_shall_cut[j+1,1]-profile_shall_cut[j,1])**2)
			slab_norm_fullstress_shall_gradient[j] = dstress_shall/dprof_shall
			# deep
			dstress_deep = (slab_norm_fullstress_deep_cut[j+1] - slab_norm_fullstress_deep_cut[j])/1.e6
			dprof_deep   = np.sqrt((profile_deep_cut[j+1,0]-profile_deep_cut[j,0])**2 + (profile_deep_cut[j+1,1]-profile_deep_cut[j,1])**2)
			slab_norm_fullstress_deep_gradient[j] = dstress_deep/dprof_deep

		elif j == (len(slab_norm_fullstress_cut)-1):
			dstress = (slab_norm_fullstress_cut[j] - slab_norm_fullstress_cut[j-1])/1.e6
			dprof   = np.sqrt((profile_cut[j,0]-profile_cut[j-1,0])**2 + (profile_cut[j,1]-profile_cut[j-1,1])**2)
			slab_norm_fullstress_gradient[j] = dstress/dprof
			# shall
			dstress_shall = (slab_norm_fullstress_shall_cut[j] - slab_norm_fullstress_shall_cut[j-1])/1.e6
			dprof_shall   = np.sqrt((profile_shall_cut[j,0]-profile_shall_cut[j-1,0])**2 + (profile_shall_cut[j,1]-profile_shall_cut[j-1,1])**2)
			slab_norm_fullstress_shall_gradient[j] = dstress_shall/dprof_shall
			# deep
			dstress_deep = (slab_norm_fullstress_deep_cut[j] - slab_norm_fullstress_deep_cut[j-1])/1.e6
			dprof_deep   = np.sqrt((profile_deep_cut[j,0]-profile_deep_cut[j-1,0])**2 + (profile_deep_cut[j,1]-profile_deep_cut[j-1,1])**2)
			slab_norm_fullstress_deep_gradient[j] = dstress_deep/dprof_deep

		else:
			dstress = (slab_norm_fullstress_cut[j+1] - slab_norm_fullstress_cut[j-1])/1.e6
			dprof   = np.sqrt((profile_cut[j+1,0]-profile_cut[j-1,0])**2 + (profile_cut[j+1,1]-profile_cut[j-1,1])**2)
			slab_norm_fullstress_gradient[j] = dstress/dprof
			# shall
			dstress_shall = (slab_norm_fullstress_shall_cut[j+1] - slab_norm_fullstress_shall_cut[j-1])/1.e6
			dprof_shall   = np.sqrt((profile_shall_cut[j+1,0]-profile_shall_cut[j-1,0])**2 + (profile_shall_cut[j+1,1]-profile_shall_cut[j-1,1])**2)
			slab_norm_fullstress_shall_gradient[j] = dstress_shall/dprof_shall
			# deep
			dstress_deep = (slab_norm_fullstress_deep_cut[j+1] - slab_norm_fullstress_deep_cut[j-1])/1.e6
			dprof_deep   = np.sqrt((profile_deep_cut[j+1,0]-profile_deep_cut[j-1,0])**2 + (profile_deep_cut[j+1,1]-profile_deep_cut[j-1,1])**2)
			slab_norm_fullstress_deep_gradient[j] = dstress_deep/dprof_deep

		# get ~central value
		if np.abs(profile_cut[j,2]-40) < dr:
			dr = np.abs(profile_cut[j,2]-40)
			slab_norm_fullstress_gradient_center = slab_norm_fullstress_gradient[j]
		if np.abs(profile_shall_cut[j,2]-40) < dr_shall:
			dr_shall = np.abs(profile_shall_cut[j,2]-40)
			slab_norm_fullstress_shall_gradient_center = slab_norm_fullstress_shall_gradient[j]
		if np.abs(profile_deep_cut[j,2]-40) < dr_deep:
			dr_deep = np.abs(profile_deep_cut[j,2]-40)
			slab_norm_fullstress_deep_gradient_center = slab_norm_fullstress_deep_gradient[j]

	print("slab-norm stress gradient (center of slab) = %.4f MPa/km" % (slab_norm_fullstress_gradient_center))

	slab_norm_fullstress_gradient_smoothed       = savgol_filter(slab_norm_fullstress_gradient,101,3)
	slab_norm_fullstress_shall_gradient_smoothed = savgol_filter(slab_norm_fullstress_shall_gradient,101,3)
	slab_norm_fullstress_deep_gradient_smoothed  = savgol_filter(slab_norm_fullstress_deep_gradient,101,3)

	# get max - min
	for j in range(len(slab_norm_fullstress_gradient_smoothed)):
		if slab_norm_fullstress_gradient_smoothed[j] > gradient_max:
			gradient_max = slab_norm_fullstress_gradient_smoothed[j]
		if slab_norm_fullstress_gradient_smoothed[j] < gradient_min:
			gradient_min = slab_norm_fullstress_gradient_smoothed[j]
	gradient_diff = gradient_max - gradient_min
	for j in range(len(slab_norm_fullstress_shall_gradient_smoothed)):
		if slab_norm_fullstress_shall_gradient_smoothed[j] > gradient_shall_max:
			gradient_shall_max = slab_norm_fullstress_shall_gradient_smoothed[j]
		if slab_norm_fullstress_shall_gradient_smoothed[j] < gradient_shall_min:
			gradient_shall_min = slab_norm_fullstress_shall_gradient_smoothed[j]
	gradient_shall_diff = gradient_shall_max - gradient_shall_min
	for j in range(len(slab_norm_fullstress_deep_gradient_smoothed)):
		if slab_norm_fullstress_deep_gradient_smoothed[j] > gradient_deep_max:
			gradient_deep_max = slab_norm_fullstress_deep_gradient_smoothed[j]
		if slab_norm_fullstress_deep_gradient_smoothed[j] < gradient_deep_min:
			gradient_deep_min = slab_norm_fullstress_deep_gradient_smoothed[j]
	gradient_deep_diff = gradient_deep_max - gradient_deep_min

	saved_stresses[ind,:] = time_dim, DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip_midmant, slab_stress_term, slab_stress_term_b, slab_stress_term_c, slabnorm_thick, Snorm_left, Snorm_right
	ind = ind + 1

	###################### plotting #########################
	fig=plt.figure()
	gs=GridSpec(5,1) 

	# plot viscosity field
	ax1=fig.add_subplot(gs[0,0])
	visc_plot = ax1.contourf(X_low/1.e3, (ymax-Y_low)/1.e3, np.log10(visc), cmap=cm.get_cmap('plasma_r'),levels=np.linspace(19,24,501))
	ymin_plot=700.e3
	ax1.set_ylim([(ymax-ymin_plot)/1.e3,0])   
	ax1.tick_params(direction='out',length=2, labelsize=6)
	ax1.annotate(''.join([str("%.1f" % (time_dim)),' Myr']), xy=(0.01,0.6), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=10,color='k')       
	ax1.annotate(''.join(['DP: mod = ',str("%.1f" % (DP_mod/1.e6)),', anal = ',str("%.1f" % (DP_anal/1.e6)),', missing =',str("%.1f" % (DP_missing/1.e6)),' MPa']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')        
	ax1.plot(llith_points[:,0], llith_points[:,1], linewidth=0.5, color='black',zorder=6)
	ax1.plot(profile_shall[:,0],profile_shall[:,1], linewidth=0.5, color='green',zorder=6)
	ax1.plot(profile_deep[:,0], profile_deep[:,1],  linewidth=0.5, color='blue',zorder=6)
	ax1.plot(profile[:,0],      profile[:,1],       linewidth=0.5, color='red',zorder=6)
	ax1.scatter(Pleft_x/1.e3,  Pleft_y/1.e3,  s=0.15,color='white',lw=0,zorder=7)
	ax1.scatter(Pright_x/1.e3, Pright_y/1.e3, s=0.15,color='white',lw=0,zorder=7)
	ax1.scatter(slab_x_left/1.e3,  slab_y_left/1.e3, s=0.15,color='yellow',lw=0,zorder=7)
	ax1.scatter(slab_x_right/1.e3, slab_y_right/1.e3, s=0.15,color='yellow',lw=0,zorder=7)
	ax1.axis('equal')

	ax3=fig.add_subplot(gs[1,0])
	ax3.set_ylabel(r"mid. $\sigma$ [MPa]",size=6.5)
	ax3.set_xlim(-10,100);
	ax3.set_ylim(-100,200); 
	ax3.tick_params(axis='x', labelsize=6)
	ax3.tick_params(axis='y', labelsize=6)
	ax3.plot(profile[:,2],sxx_slab/1.e6,label='sxx',color='tan', linestyle='-', linewidth=1.5, zorder=2)
	ax3.plot(profile[:,2],sxy_slab/1.e6,label='sxy',color='indianred',      linestyle='-', linewidth=1.5, zorder=2)
	ax3.plot(profile[:,2],-1.*sxx_slab/1.e6,label='syy',color='plum', linestyle='-', linewidth=1.5, zorder=2)
	ax3.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)

	# plot shear stress in slab
	ax4=fig.add_subplot(gs[2,0])
	ax4.set_ylabel(r"$\tau_{n}$ [MPa]",size=6.5)
	ax4.set_xlim(-10,100); 
	ax4.set_ylim(-25,10); 
	ax4.tick_params(axis='x', labelsize=6)
	ax4.tick_params(axis='y', labelsize=6)
	ax4.plot(profile[:,2],      slab_shear_stress/1.e6,      color='red', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
	ax4.plot(profile_shall[:,2],slab_shear_stress_shall/1.e6, color='green', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
	ax4.plot(profile_deep[:,2], slab_shear_stress_deep/1.e6,  color='blue', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
	ax4.plot(profile_cut[:,2],      slab_shear_stress_cut/1.e6,       color='red', linestyle='-', linewidth=1.5,zorder=3)  
	ax4.plot(profile_shall_cut[:,2],slab_shear_stress_shall_cut/1.e6, color='green', linestyle='-', linewidth=1.5,zorder=3)  
	ax4.plot(profile_deep_cut[:,2], slab_shear_stress_deep_cut/1.e6,  color='blue', linestyle='-', linewidth=1.5,zorder=3) 
	ax4.annotate(''.join(['slab stress term = ',str("%.1f" % (-1.0*slab_stress_term/1.e6)),', dev stress term = ',str("%.1f" % (Snorm_contrib/1.e6)),' misfit = ',str("%.1f" % (((Snorm_contrib-slab_stress_term)-DP_missing)/1.e6)),' MPa']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')        
	ax4.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)

	# plot norm stress (in slab)
	ax5=fig.add_subplot(gs[3,0])
	ax5.set_ylabel(r"$\sigma_{n}$ [MPa]",size=6.5)
	ax5.set_xlim(-10,100); 
	ax5.set_ylim(-45,10); 
	ax5.tick_params(axis='x', labelsize=6)
	ax5.tick_params(axis='y', labelsize=6)
	ax5.plot(profile[:,2],  slab_norm_fullstress/1.e6,  color='red', linestyle='-', linewidth=0.5,zorder=2,alpha=0.5)  
	ax5.plot(profile_cut[:,2],  slab_norm_fullstress_cut/1.e6,  color='red', linestyle='-', linewidth=1,zorder=3)  
	ax5.plot(profile_shall[:,2],  slab_norm_fullstress_shall/1.e6,  color='green', linestyle='-', linewidth=0.5,zorder=2,alpha=0.5)  
	ax5.plot(profile_shall_cut[:,2],  slab_norm_fullstress_shall_cut/1.e6,  color='green', linestyle='-', linewidth=1,zorder=3) 
	ax5.plot(profile_deep[:,2],  slab_norm_fullstress_deep/1.e6,  color='blue', linestyle='-', linewidth=0.5,zorder=2,alpha=0.5)  
	ax5.plot(profile_deep_cut[:,2],  slab_norm_fullstress_deep_cut/1.e6,  color='blue', linestyle='-', linewidth=1,zorder=3)  

	# plot norm stress gradient (in slab)
	ax6=fig.add_subplot(gs[4,0])
	ax6.set_ylabel(r"$\sigma_{n}$/dr [MPa/km]",size=6.5)
	ax6.set_xlim(-10,100); 
	ax6.set_ylim(-0.6,0.6); 
	ax6.tick_params(axis='x', labelsize=6)
	ax6.tick_params(axis='y', labelsize=6)
	ax6.plot(profile_cut[:,2],       slab_norm_fullstress_gradient-slab_norm_fullstress_gradient_center,  color='red', linestyle='-', linewidth=0.5,alpha=0.5,zorder=2)  
	ax6.plot(profile_cut[:,2],       slab_norm_fullstress_gradient_smoothed-slab_norm_fullstress_gradient_center,  color='red', linestyle='-', linewidth=1,zorder=3)
	ax6.plot(profile_shall_cut[:,2],     slab_norm_fullstress_shall_gradient-slab_norm_fullstress_shall_gradient_center,  color='green', linestyle='-', linewidth=0.5,alpha=0.5,zorder=2)  
	ax6.plot(profile_shall_cut[:,2], slab_norm_fullstress_shall_gradient_smoothed-slab_norm_fullstress_shall_gradient_center,  color='green', linestyle='-', linewidth=1,zorder=3)  
	ax6.plot(profile_deep_cut[:,2],      slab_norm_fullstress_deep_gradient-slab_norm_fullstress_deep_gradient_center,  color='blue', linestyle='-', linewidth=0.5,alpha=0.5,zorder=2)  
	ax6.plot(profile_deep_cut[:,2],  slab_norm_fullstress_deep_gradient_smoothed-slab_norm_fullstress_deep_gradient_center,  color='blue', linestyle='-', linewidth=1,zorder=3)    
	ax6.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)
	ax6.annotate(''.join(['gradient diffs, center = ',str("%.2f" % (gradient_diff)),', shall = ',str("%.2f" % (gradient_shall_diff)),', deep = ',str("%.2f" % (gradient_deep_diff)),' MPa/km']), xy=(0.01,0.8), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')        
	normalize_factor       = (slabnorm_thick/1.e3)/(DP_mod/1.e6)# km/MPa 
	ax6.annotate(''.join(['normalized, center = ',str("%.2f" % (gradient_diff/normalize_factor)),', shall = ',str("%.2f" % (gradient_shall_diff/normalize_factor)),', deep = ',str("%.2f" % (gradient_deep_diff/normalize_factor)),' MPa/km']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')        


	# print("saving fields figure to %s..." % plotname)
	plt.savefig(plotname, bbox_inches='tight', format='png', dpi=1000)
	plt.clf()


##############################################################################################################################
# simple scatter plot
plot_name_png=''.join([plot_loc,'/z',str(analysis_depth/1.e3),'km.scatter.png'])

fig=plt.figure()
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

# plot viscosity field
ax=fig.add_subplot(gs[0,0])

plt.scatter(saved_stresses[:,4]/1.e6,saved_stresses[:,3]/1.e6,color='blue',s=8,lw=0.0,zorder=3,alpha=0.6) # uncorrected
plt.scatter(saved_stresses[:,4]/1.e6,(saved_stresses[:,3]-saved_stresses[:,6])/1.e6,color='blue',s=10,edgecolor='black',label='230 km',lw=0.25,zorder=4) # corrected
for i in range(len(saved_stresses)):
	plt.plot([saved_stresses[i,4]/1.e6, saved_stresses[i,4]/1.e6], [saved_stresses[i,3]/1.e6, (saved_stresses[i,3]-saved_stresses[i,6])/1.e6], color='blue', linewidth=1, zorder=2, alpha=0.5)

# axis stuff
plt.xlim(-17.5,  40); plt.ylim(-17.5,  40)
plt.plot([-17.5, 40], [-17.5, 40], color='black', linewidth=1, zorder=1)
plt.axhline(y=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
plt.axvline(x=0, color='gray',linestyle='--',linewidth=0.5, zorder=0)
ax.tick_params(axis='x', labelsize=6)
ax.tick_params(axis='y', labelsize=6)
plt.xlabel("slab buoyancy  [MPa]",size=7.5)
plt.ylabel("extracted stress  [MPa]",size=7.5)
fixed_aspect_ratio(1)


plt.savefig(plot_name_png, bbox_inches='tight', format='png', dpi=500)


