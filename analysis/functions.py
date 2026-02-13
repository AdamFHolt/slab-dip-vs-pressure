#!/bin/python
import numpy as np
import scipy.special
import scipy.integrate
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import sys, os, math, statistics

def create_grid(xmin,xmax,ymin,ymax,grid_res):

	x_low = np.linspace(xmin,xmax,int((xmax-xmin)/grid_res))
	y_low =  np.linspace(ymin,ymax,int((ymax-ymin)/grid_res))
	X_grid, Y_grid = np.meshgrid(x_low,y_low)	

	return X_grid, Y_grid


def get_slab_midplane(llith_points_tmp,cutoff_shall,cutoff_deep):
	
	# remove shallow and deep points
	n = 0
	for i in range(len(llith_points_tmp)):
		if llith_points_tmp[i,1] >= cutoff_shall and llith_points_tmp[i,1] <= cutoff_deep:
			n = n + 1
	llith_points_tmp2 = np.zeros([n,2])
	ind = 0
	for i in range(len(llith_points_tmp)):
		if llith_points_tmp[i,1] >= cutoff_shall and llith_points_tmp[i,1] <= cutoff_deep:
				llith_points_tmp2[ind,:] = llith_points_tmp[i,:]
				ind = ind + 1
	llith_points_tmp = llith_points_tmp2;
	# remove the other plane (i.e., base of slab)
	depth_old = 0.
	for j in range(len(llith_points_tmp)):
		depth = llith_points_tmp[j,1] 
		if depth < depth_old:
			ind = j
			break 
		depth_old = depth
	llith_points = np.delete(llith_points_tmp, range(0,ind,1), axis=0)
	llith_points = np.flipud(llith_points)
	# add along-slab distance
	along_slab_dist = np.zeros((len(llith_points),1))
	for i in range(1,len(llith_points)):
		along_slab_dist[i,0] = along_slab_dist[i-1,0] + np.sqrt((llith_points[i,0]-llith_points[i-1,0])**2 + (llith_points[i,1]-llith_points[i-1,1])**2)
	llith_points = np.concatenate((llith_points,along_slab_dist), axis=1)

	return llith_points


def extract_horiz_prof(model_data,depth,ymax,dz,x_col,y_col):

	# get horiz profile in mid-mantle 
	midmant_prof_loc = ymax - depth                          
	midmant_prof = model_data[model_data[:,y_col] < (midmant_prof_loc+dz)] 
	midmant_prof = midmant_prof[midmant_prof[:,y_col] > (midmant_prof_loc-dz)] 
	midmant_prof = midmant_prof[midmant_prof[:,x_col].argsort()] # sort by x

	return midmant_prof

def get_slablocation_from_horiz_prof(midmant_prof,x_col,y_col,c_ulith_col,c_llith_col,ymax):

	# get slab location
	x_left = 1e9; x_right = 0; x_center = 0;
	for i in range(len(midmant_prof)):
		x = midmant_prof[i,x_col]
		y = (ymax - midmant_prof[i,y_col])/1.e3
		c_ulith = midmant_prof[i,c_ulith_col]
		c_llith = midmant_prof[i,c_llith_col]
		if x < x_left and c_llith > 0.5:
			x_left = x
			y_left = y
			ind_left = i
		if x > x_center and c_llith > 0.5:
			x_center = x
			y_center = y
			ind_center = i 
		if x > x_right and c_ulith > 0.5:
			x_right = x
			y_right = y
			ind_right = i

	return x_left, y_left, x_center, y_center, x_right, y_right


def get_slabvisc_from_horiz_prof(midmant_prof,x_col,y_col,visc_col,x_center,y_center,ymax):

	# get slab location
	dist_to_cent = 1e12
	for i in range(len(midmant_prof)):
		x = midmant_prof[i,x_col]            # m
		y = ymax - midmant_prof[i,y_col]     # m
		dist_to_cent_tmp = np.sqrt((x-x_center)**2 + (y-(y_center*1.e3))**2)
		if dist_to_cent_tmp < dist_to_cent:
			dist_to_cent = dist_to_cent_tmp
			slab_viscosity   = midmant_prof[i,visc_col]
			slab_viscosity_x = x
			slab_viscosity_y = y
	return slab_viscosity, slab_viscosity_x/1e3, slab_viscosity_y/1e3


def get_platevels_from_horiz_prof(surf_prof,x_col,y_col,c_crust_col,vx_col,ymax):
   
	# get trench location
	trench_loc = 0;
	for i in range(len(surf_prof)):
		if (surf_prof[i,c_crust_col] > 0.5 and surf_prof[i,c_crust_col] > trench_loc):
			trench_loc = surf_prof[i,x_col]

	# get plate velocities either side of trench
	vsp_tot = 0; nsp = 0
	vop_tot = 0; nop = 0
	for i in range(len(surf_prof)):
		if surf_prof[i,x_col] > (trench_loc - 250e3) and surf_prof[i,x_col] < (trench_loc - 200e3):
			vsp_tot = vsp_tot + surf_prof[i,vx_col]
			nsp = nsp + 1
		if surf_prof[i,x_col] > (trench_loc + 200e3) and surf_prof[i,x_col] < (trench_loc + 250e3):
			vop_tot = vop_tot + surf_prof[i,vx_col]
			nop = nop + 1
	vsp = 100.*(vsp_tot/nsp)
	vop = 100.*(vop_tot/nop)
	vc  = vsp - vop

	return vc, vsp, vop, trench_loc


def get_slablocation_from_horiz_prof_coremodels(midmant_prof,x_col,y_col,c_lith_col,c_lcore_col,ymax):

	# get slab location
	x_left = 1e9; x_right = 0; x_center = 0;
	for i in range(len(midmant_prof)):
		x = midmant_prof[i,x_col]
		y = (ymax - midmant_prof[i,y_col])/1.e3
		c_lith  = midmant_prof[i,c_lith_col]
		c_lcore = midmant_prof[i,c_lcore_col]
		if x < x_left and c_lith > 0.5:
			x_left = x
			y_left = y
			ind_left = i
		if x > x_center and c_lcore > 0.5:
			x_center = x
			y_center = y
			ind_center = i 
		if x > x_right and c_lith > 0.5:
			x_right = x
			y_right = y
			ind_right = i

	return x_left, y_left, x_center, y_center, x_right, y_right

def get_dip_slab_midplane(llith_points):

	# get dip of whole surface 
	dips = np.zeros((len(llith_points),1))
	for i in range(len(llith_points)):
		
		if i == 0:
			dx = llith_points[i+1,0] - llith_points[i,0]
			dy = llith_points[i+1,1] - llith_points[i,1]   
		elif i == len(llith_points)-1:
			dx = llith_points[i,0] - llith_points[i-1,0]
			dy = llith_points[i,1] - llith_points[i-1,1]
		else:
			dx = llith_points[i+1,0]   - llith_points[i-1,0]
			dy = llith_points[i+1,1]   - llith_points[i-1,1] 

		if dx < 0:
			dx = -1. * dx
			dips[i] = 180. - np.rad2deg(np.arctan(dy/dx))
		else:
			dips[i] = np.rad2deg(np.arctan(dy/dx))

	return dips 


def get_curvature_slab_midplane(llith_points_tmp,dips):

	# get curvature 
	K_unsmoothed = np.zeros((len(llith_points_tmp),1))

	for i in range(len(K_unsmoothed)):
		
		if i == 0:
			ds   = llith_points_tmp[i+1,2] - llith_points_tmp[i,2]  # km
			ddip = dips[i+1] - dips[i]                              # degs
		elif i == len(llith_points_tmp)-1:
			ds   = llith_points_tmp[i,2] - llith_points_tmp[i-1,2]
		else:
			ds   = llith_points_tmp[i+1,2]   - llith_points_tmp[i-1,2]
			ddip = dips[i+1]   - dips[i-1] 

		K_unsmoothed[i] = np.deg2rad(ddip)/(ds*1.e3)                # rad/m

	K = savgol_filter(K_unsmoothed[:,0],601,3)

	# get dcurvature/ds 
	dKds_unsmoothed = np.zeros((len(llith_points_tmp),1))
	for i in range(len(dKds_unsmoothed)):
		
		if i == 0:
			ds = llith_points_tmp[i+1,2] - llith_points_tmp[i,2]   # km
			dK = K[i+1] - K[i]                                     # rad/m
		elif i == len(llith_points_tmp)-1:
			ds = llith_points_tmp[i,2] - llith_points_tmp[i-1,2]
			dK = K[i] - K[i-1]
		else:
			ds = llith_points_tmp[i+1,2]   - llith_points_tmp[i-1,2]
			dK = K[i+1]   - K[i-1] 

		dKds_unsmoothed[i] = dK/(ds*1.e3)                          # rad/m^2

	dKds = savgol_filter(dKds_unsmoothed[:,0],601,3)

	return K, dKds, K_unsmoothed, dKds_unsmoothed


def get_dip_at_certain_depth(dips,llith_points,analysis_depth):

	dz = 1e9
	for i in range(1,len(llith_points)):
		if np.abs(llith_points[i,1]-(analysis_depth/1.e3)) < dz:
			dz = np.abs(llith_points[i,1]-(analysis_depth/1.e3))
			dip=dips[i]
			xloc=llith_points[i,0]
			sloc=llith_points[i,2]

	return dip, xloc, sloc


def get_stress_at_certain_depth(bending_stress,llith_points,analysis_depth):

	dz = 1e9
	for i in range(1,len(llith_points)):
		if np.abs(llith_points[i,1]-(analysis_depth/1.e3)) < dz:
			dz = np.abs(llith_points[i,1]-(analysis_depth/1.e3))
			stress=bending_stress[i]
			xloc=llith_points[i,0]
			sloc=llith_points[i,2]

	return stress, xloc, sloc


def get_curvature_at_certain_depth(K,dK,llith_points,analysis_depth):

	dz = 1e9
	for i in range(1,len(llith_points)):
		if np.abs(llith_points[i,1]-(analysis_depth/1.e3)) < dz:
			dz = np.abs(llith_points[i,1]-(analysis_depth/1.e3))
			K_midmant=K[i]
			dK_midmant=dK[i]
			xloc=llith_points[i,0]
			sloc=llith_points[i,2]

	return K_midmant, dK_midmant, xloc, sloc


def get_stress_profile(model_data,slab_x_left,slab_y_left,slab_x_right,slab_y_right,ymax,x_col,y_col,sxx_col,sxy_col,P_col):

	profile = np.zeros((200,3))
	for j in range(len(profile)):
		x = (slab_x_left/1.e3) + (j/(len(profile)-1)) * (slab_x_right-slab_x_left)/1.e3 
		y = (slab_y_left/1.e3) + (j/(len(profile)-1)) * (slab_y_right-slab_y_left)/1.e3
		dist = np.sqrt((x-(slab_x_left/1.e3))**2 + (y-(slab_y_left/1.e3))**2)
		profile[j,:] = x, y, dist
	sxx_slab = griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,sxx_col], (profile[:,0], profile[:,1]), method='linear')
	sxy_slab = griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,sxy_col], (profile[:,0], profile[:,1]), method='linear')
	P_slab   = griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,P_col], (profile[:,0], profile[:,1]), method='linear')

	return sxx_slab, sxy_slab, P_slab, profile


def convert_to_slabnorm_shearstress(sxx_slab,sxy_slab,dip_midmant):

	slab_shear_stress = np.zeros((len(sxx_slab),1))
	for i in range(0,len(sxx_slab)):
		norm_contrib = 2.*sxx_slab[i]*np.sin(np.deg2rad(dip_midmant))*np.cos(np.deg2rad(dip_midmant));
		shear_contrib = sxy_slab[i]*((np.cos(np.deg2rad(dip_midmant)))**2 - (np.sin(np.deg2rad(dip_midmant)))**2);
		slab_shear_stress[i] = norm_contrib + shear_contrib

	return slab_shear_stress

# def convert_to_slabnorm_normstress(sxx_slab,sxy_slab, P_slab, dip_midmant):

# 	slab_norm_devstress  = np.zeros((len(sxx_slab),1))
# 	slab_norm_fullstress = np.zeros((len(sxx_slab),1))
# 	for i in range(0,len(sxx_slab)):
# 		norm_contrib        = sxx_slab[i]*((np.sin(np.deg2rad(dip_midmant)))**2 - (np.cos(np.deg2rad(dip_midmant)))**2)
# 		shear_contrib       = sxy_slab[i]*np.sin(np.deg2rad(2.*dip_midmant))
# 		slab_norm_devstress[i]  = norm_contrib + shear_contrib
# 		slab_norm_fullstress[i] = slab_norm_devstress[i] + P_slab[i]

# 	return slab_norm_devstress, slab_norm_fullstress


def convert_to_slabnorm_normstress(sxx_slab,sxy_slab, P_slab, dip_midmant):

	slab_norm_devstress  = np.zeros((len(sxx_slab),1))
	slab_norm_fullstress = np.zeros((len(sxx_slab),1))
	for i in range(0,len(sxx_slab)):
		norm_contrib        = sxx_slab[i]*((np.cos(np.deg2rad(dip_midmant)))**2 - (np.sin(np.deg2rad(dip_midmant)))**2)
		shear_contrib       = -1.0 * sxy_slab[i]*np.sin(np.deg2rad(2.*dip_midmant))
		slab_norm_devstress[i]  = norm_contrib + shear_contrib
		slab_norm_fullstress[i] = slab_norm_devstress[i] + P_slab[i]

	return slab_norm_devstress, slab_norm_fullstress

def convert_to_slabnorm_normstress_TEST(sxx_slab,sxy_slab, P_slab, dip_midmant):

	slab_norm_devstress_test  = np.zeros((len(sxx_slab),1))
	slab_norm_fullstress_test = np.zeros((len(sxx_slab),1))
	for i in range(0,len(sxx_slab)):
		norm_contrib        = sxx_slab[i]*((np.cos(np.deg2rad(dip_midmant)))**2 - (np.sin(np.deg2rad(dip_midmant)))**2)
		shear_contrib       = sxy_slab[i]*np.sin(np.deg2rad(2.*dip_midmant))  # NOTE: this has lost the minus sign
		slab_norm_devstress_test[i]  = norm_contrib + shear_contrib
		slab_norm_fullstress_test[i] = slab_norm_devstress_test[i] + P_slab[i]

	return slab_norm_devstress_test, slab_norm_fullstress_test	


def get_farfieldP(y_center,ymax,dz,model_data,P_col,x_col,y_col):

    prof_loc = ymax - (y_center*1.e3)
    prof = model_data[model_data[:,y_col] < (prof_loc+dz)]
    prof = prof[prof[:,y_col] > (prof_loc-dz)]
    prof = prof[prof[:,x_col].argsort()] # sort by x

    min_loc = 1e9
    for i in range(len(prof)):
        x = prof[i,x_col]
        if x < min_loc:
            min_loc = x
            Pref = prof[i,P_col]

    return Pref



def get_nearslab_stresses(y_center,slabnorm_thick_lower,slabnorm_thick_upper,dip_midmant,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz):
	z_shallow = y_center - (slabnorm_thick_lower*1.e-3)*np.cos(np.deg2rad(dip_midmant)) # km
	z_deep    = y_center + (slabnorm_thick_upper*1.e-3)*np.cos(np.deg2rad(dip_midmant)) # km
	shallow_prof_loc = ymax - (z_shallow*1.e3)
	shallow_prof = model_data[model_data[:,y_col] < (shallow_prof_loc+dz)] 
	shallow_prof = shallow_prof[shallow_prof[:,y_col] > (shallow_prof_loc-dz)] 
	shallow_prof = shallow_prof[shallow_prof[:,x_col].argsort()] # sort by x
	deep_prof_loc = ymax - (z_deep*1.e3)                          
	deep_prof = model_data[model_data[:,y_col] < (deep_prof_loc+dz)] 
	deep_prof = deep_prof[deep_prof[:,y_col] > (deep_prof_loc-dz)] 
	deep_prof = deep_prof[deep_prof[:,x_col].argsort()] # sort by x

	# get point locations for dp
	slab_x_right = 0;
	for i in range(len(shallow_prof)):
		x = shallow_prof[i,x_col]
		y = (ymax - shallow_prof[i,y_col])
		c_ulith = shallow_prof[i,c_ulith_col]
		if x > slab_x_right and c_ulith > 0.5:
			slab_x_right = x
			slab_y_right = y
			ind_right = i
	slab_x_left = 1e9;
	for i in range(len(deep_prof)):
		x = deep_prof[i,x_col]
		y = (ymax - deep_prof[i,y_col])
		c_llith = deep_prof[i,c_llith_col]
		if x < slab_x_left and c_llith > 0.5:
			slab_x_left = x
			slab_y_left = y
			ind_left = i

	# find pressure at these points
	misfit_left = 1.e9
	misfit_right = 1.e9
	pos_left = slab_x_left - ds
	pos_right = slab_x_right + ds
	for i in range(len(shallow_prof)):
		x = shallow_prof[i,x_col]
		y = ymax - shallow_prof[i,y_col]
		misfit_right_tmp = np.abs(x-pos_right)
		if misfit_right_tmp < misfit_right:
			misfit_right = misfit_right_tmp
			Pright       = shallow_prof[i,P_col] 
			Sxxright     = shallow_prof[i,sxx_col]
			Sxyright     = shallow_prof[i,sxy_col]
			Pright_x     = shallow_prof[i,x_col]
			Pright_y	 = ymax - shallow_prof[i,y_col]

	for i in range(len(deep_prof)):
		x = deep_prof[i,x_col]
		y = ymax - deep_prof[i,y_col]
		misfit_left_tmp  = np.abs(x-pos_left)
		if misfit_left_tmp < misfit_left:
			misfit_left = misfit_left_tmp
			Pleft       = deep_prof[i,P_col] 
			Sxxleft     = deep_prof[i,sxx_col] 
			Sxyleft     = deep_prof[i,sxy_col] 
			Pleft_x     = deep_prof[i,x_col]
			Pleft_y		= ymax - deep_prof[i,y_col]

	return Pleft, Pright, Sxxleft, Sxxright, Sxyleft, Sxyright, Pleft_x, Pright_x, Pleft_y, Pright_y, slab_x_left, slab_y_left, slab_x_right, slab_y_right



def get_average_stress_and_dips(model_file,drho):

	tot_bendingstress = 0; 
	tot_dip = 0; 
	tot_DP  = 0
	tot_DPanal1 = 0; 
	tot_DP2 = 0; 
	for i in range(2,len(model_file)):

		# bending stresses
		tot_bendingstress = tot_bendingstress + (model_file[i,6]/1.e6) # MPa
		# dips
		tot_dip = tot_dip + model_file[i,5] 	 # degrees
		# DP
		tot_DP = tot_DP + (model_file[i,3]/1.e6) # MPa
		# DP2 (with bending stress!)
		tot_DP2 = tot_DP2 + (model_file[i,3]/1.e6) - (model_file[i,6]/1.e6)
		# DP from dip
		tot_DPanal1 = tot_DPanal1 + ((drho * 9.81 * model_file[i,9] * np.cos(np.deg2rad(model_file[i,5])))/1.e6) # MPa

	num = len(model_file) - 2
	mean_bending 	= tot_bendingstress/num
	mean_dip 		= tot_dip/num
	mean_DP 		= tot_DP/num
	mean_DPanal1 	= tot_DPanal1/num
	mean_DP2 		= tot_DP2/num

	return mean_bending,mean_dip,mean_DP,mean_DPanal1,mean_DP2


# def get_stress_gradient(stress_profile,profile):

# 	# normal stress gradient calculation:
# 	stress_profile_gradient = np.zeros(len(stress_profile))
# 	for j in range(len(stress_profile)):
# 		if j == 0:
# 			dstress = (stress_profile[j+1] - stress_profile[j])/1.e6
# 			dprof   = np.sqrt((profile[j+1,0]-profile[j,0])**2 + (profile[j+1,1]-profile[j,1])**2)
# 			stress_profile_gradient[j] = dstress/dprof  # MPa / km
# 		elif j == (len(stress_profile)-1):
# 			dstress = (stress_profile[j] - stress_profile[j-1])/1.e6
# 			dprof   = np.sqrt((profile[j,0]-profile[j-1,0])**2 + (profile[j,1]-profile[j-1,1])**2)
# 			stress_profile_gradient[j] = dstress/dprof
# 		else:
# 			dstress = (stress_profile[j+1] - stress_profile[j-1])/1.e6
# 			dprof   = np.sqrt((profile[j+1,0]-profile[j-1,0])**2 + (profile[j+1,1]-profile[j-1,1])**2)
# 			stress_profile_gradient[j] = dstress/dprof

# 	stress_profile_gradient_smoothed = savgol_filter(stress_profile_gradient,101,3)
# Play all]
# 		if stress_profile_gradient_smoothed[j] < gradient_min:
# 			gradient_min = stress_profile_gradient_smoothed[j]
# 	gradient_diff = gradient_max - gradient_min

# 	return stress_profile_gradient_smoothed, gradient_diff


def get_misfit_mean_and_stdev(mod,tmin):

	mod_analyze = mod[tmin:,:]
	misfit_total_woshear = 0; misfit_total_wshear = 0; n = 0
	array_for_stdev = np.zeros((len(mod_analyze),2))
	for i in range(len(mod_analyze)):
		misfit_total_woshear = misfit_total_woshear + (mod_analyze[i,3]-mod_analyze[i,4])/1.e6
		misfit_total_wshear  = misfit_total_wshear  + (mod_analyze[i,3]-mod_analyze[i,4]-mod_analyze[i,6]-mod_analyze[i,17])/1.e6
		array_for_stdev[i,0] = (mod_analyze[i,3]-mod_analyze[i,4])/1.e6
		array_for_stdev[i,1] = (mod_analyze[i,3]-mod_analyze[i,4]-mod_analyze[i,6]-mod_analyze[i,17])/1.e6
		n = n + 1
	misfit_mean_woshear  = misfit_total_woshear/n 
	misfit_mean_wshear   = misfit_total_wshear/n
	misfit_stdev_woshear = statistics.stdev(array_for_stdev[:,0])
	misfit_stdev_wshear  = statistics.stdev(array_for_stdev[:,1])

	return np.abs(misfit_mean_woshear), np.abs(misfit_mean_wshear), misfit_stdev_woshear, misfit_stdev_wshear


def get_misfit_mean_and_stdev_nondim(mod,tmin):

	mod_analyze = mod[tmin:,:]
	misfit_total_woshear = 0; misfit_total_wshear = 0; n = 0
	array_for_stdev = np.zeros((len(mod_analyze),2))
	for i in range(len(mod_analyze)):
		misfit_total_woshear = misfit_total_woshear + 100.*np.abs((mod_analyze[i,3]-mod_analyze[i,4])/mod_analyze[i,4])
		misfit_total_wshear  = misfit_total_wshear  + 100.*np.abs((mod_analyze[i,3]-mod_analyze[i,4]-mod_analyze[i,6]-mod_analyze[i,17])/mod_analyze[i,4])
		array_for_stdev[i,0] = 100.*(mod_analyze[i,3]-mod_analyze[i,4])/mod_analyze[i,4]
		array_for_stdev[i,1] = 100.*(mod_analyze[i,3]-mod_analyze[i,4]-mod_analyze[i,6]-mod_analyze[i,17])/mod_analyze[i,4]
		n = n + 1
	misfit_mean_woshear  = misfit_total_woshear/n 
	misfit_mean_wshear   = misfit_total_wshear/n
	misfit_stdev_woshear = statistics.stdev(array_for_stdev[:,0])
	misfit_stdev_wshear  = statistics.stdev(array_for_stdev[:,1])

	return misfit_mean_woshear, misfit_mean_wshear, misfit_stdev_woshear, misfit_stdev_wshear


def get_avg_forces_nondim(mod,tmin):

	mod_analyze = mod[tmin:,:]
	total_anal = 0; total_mod = 0
	total_DP    = 0
	total_shear = 0
	total_norm  = 0
	array_stdev   = np.zeros((len(mod_analyze),1))
	n = 0
	for i in range(len(mod_analyze)):
		total_anal        = total_anal + mod_analyze[i,4]
		mod_force		  = mod_analyze[i,3]-mod_analyze[i,6]+mod_analyze[i,17]
		total_mod         = total_mod  + mod_force
		total_DP    = total_DP + mod_analyze[i,3]
		total_shear = total_shear + ((-1.0 * mod_analyze[i,6])) 
		total_norm  = total_norm  + ((mod_analyze[i,17]))
		array_stdev[i,0] = mod_force/mod_analyze[i,4]
		n = n + 1

	avg_anal 		= total_anal/n
	avg_mod  		= total_mod/n
	avg_DP	    	= total_DP/n
	avg_shear 		= total_shear/n
	avg_norm  		= total_norm/n

	return avg_anal, avg_mod, avg_DP, avg_shear, avg_norm, statistics.stdev(array_stdev[:,0])

def get_avg_forces_nondim_curvethresh(mod,tmin,curve_thresh):

	mod_analyze = mod[tmin:,:]
	total_anal = 0; total_mod = 0
	total_DP    = 0
	total_shear = 0
	total_norm  = 0
	total_K_thresh = 0
	total_K_full = 0
	array_stdev   = np.zeros((len(mod_analyze),1))
	n = 0; nfull = 0
	for i in range(len(mod_analyze)):
		total_K_full = total_K_full + (mod_analyze[i,11]*1e3)
		nfull = nfull + 1
		if (mod_analyze[i,11]*1e3) < curve_thresh:
			total_anal        = total_anal + mod_analyze[i,4]
			mod_force		  = mod_analyze[i,3]-mod_analyze[i,6]+mod_analyze[i,17]
			total_mod         = total_mod  + mod_force
			total_DP    = total_DP + mod_analyze[i,3]
			total_shear = total_shear + ((-1.0 * mod_analyze[i,6])) 
			total_norm  = total_norm  + ((mod_analyze[i,17]))
			total_K_thresh = total_K_thresh + (mod_analyze[i,11]*1e3)
			array_stdev[i,0] = mod_force/mod_analyze[i,4]
			n = n + 1

	K_full = total_K_full/nfull

	if n == 0:
		# print(mod)
		return 0, 0, 0, 0, 0, 0, 0, K_full
	else:
		avg_anal 		= total_anal/n
		avg_mod  		= total_mod/n
		avg_DP	    	= total_DP/n
		avg_shear 		= total_shear/n
		avg_norm  		= total_norm/n
		avg_K_thresh 	= total_K_thresh/n

	return avg_anal, avg_mod, avg_DP, avg_shear, avg_norm, statistics.stdev(array_stdev[:,0]), avg_K_thresh, K_full



def get_avg_forces_curvethresh(mod,tmin,curve_thresh):

	mod_analyze = mod[tmin:,:]
	total_anal = 0; total_mod = 0
	total_DP    = 0
	total_shear = 0
	total_norm  = 0
	total_K_thresh = 0
	total_K_full = 0
	array_stdev   = np.zeros((len(mod_analyze),1))
	n = 0; nfull = 0
	for i in range(len(mod_analyze)):
		total_K_full = total_K_full + (mod_analyze[i,11]*1e3)
		nfull = nfull + 1
		if (mod_analyze[i,11]*1e3) < curve_thresh:
			total_anal        = total_anal + mod_analyze[i,4]
			mod_force		  = mod_analyze[i,3]-mod_analyze[i,6]+mod_analyze[i,17]
			total_mod         = total_mod  + mod_force
			total_DP    = total_DP + mod_analyze[i,3]
			total_shear = total_shear + ((-1.0 * mod_analyze[i,6])) 
			total_norm  = total_norm  + ((mod_analyze[i,17]))
			total_K_thresh = total_K_thresh + (mod_analyze[i,11]*1e3)
			array_stdev[i,0] = mod_force/mod_analyze[i,4]
			n = n + 1

	K_full = total_K_full/nfull

	if n == 0:
		# print(mod)
		return 0, 0, 0, 0, 0, 0, 0, K_full
	else:
		avg_anal 		= total_anal/n
		avg_mod  		= total_mod/n
		avg_DP	    	= total_DP/n
		avg_shear 		= total_shear/n
		avg_norm  		= total_norm/n
		avg_K_thresh 	= total_K_thresh/n

	return avg_anal, avg_mod, avg_DP, avg_shear, avg_norm, statistics.stdev(array_stdev[:,0]), avg_K_thresh, K_full


def get_curvature_mean_and_stdev(mod,tmin):

	mod_analyze = mod[tmin:,:]
	total = 0; n = 0
	array_for_stdev = np.zeros((len(mod_analyze),1))
	for i in range(len(mod_analyze)):
		total  = total + (mod_analyze[i,11]*1000.)
		array_for_stdev[i,0] = mod_analyze[i,11]*1000.
		n = n + 1
	mean  = total/n
	stdev = statistics.stdev(array_for_stdev[:,0])

	return np.abs(mean), stdev
	
