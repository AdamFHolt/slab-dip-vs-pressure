#!/bin/python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator
import sys, os, subprocess
from scipy.signal import savgol_filter
from scipy.interpolate import splrep, splev, UnivariateSpline

model_name=str(sys.argv[1])	            
max_time=int(sys.argv[2])   # largest number in csv_outputs/ filenames
analysis_depth = 330.e3 # km

xmax=5800.e3
ymax=1450.e3
mant_visc = 2.5e20
slab_visc = mant_visc * 500.
drho = 50. # kg/m3

# ASPECT output 
csvs_loc =  'csv_outputs/'
models_loc =  'raw_outputs/'
stats_file = ''.join([models_loc,str(model_name),'/statistics'])
model_output_dt  = 50 # output dt as set in ASPECT .prm file (for getting the dimensional time)
num_header_lines = 16 # num header lines in stats_files (for getting the dimensional time)

# where to put the plots
plot_loc   = ''.join(['plots/evolution/',str(model_name)])
if not os.path.exists(plot_loc):
	os.mkdir(plot_loc)

# column numbers of the relevant properties in the .csv file. 
visc_col=26; 	vx_col=0;  		 	vy_col=1;
c_crust_col=23; c_ulith_col = 24; 	c_llith_col = 25;
P_col = 29; 	x_col=30;  		 	y_col = 31;
sxx_col = 3;	syy_col = 7;		sxy_col = 4; 

# create low res grid to interpolate stuff onto (for plotting)
xmin_plot=0; ymin_plot=500.e3; grid_res=2.5e3
x_low = np.linspace(xmin_plot,xmax,int((xmax-xmin_plot)/grid_res))
y_low =  np.linspace(ymin_plot,ymax,int((ymax-ymin_plot)/grid_res))
X_low, Y_low = np.meshgrid(x_low,y_low)	
# higher res grid
xmin_plot2=2500.e3; xmax_plot2=3750.e3
ymin_plot2=ymax-600.e3; grid_res2=0.5e3
x_low2 = np.linspace(xmin_plot2,xmax_plot2,int((xmax_plot2-xmin_plot2)/grid_res2))
y_low2 =  np.linspace(ymin_plot2,ymax,int((ymax-ymin_plot2)/grid_res2))
X_low2, Y_low2 = np.meshgrid(x_low2,y_low2)	

# for time in range(8,max_time+1,1):
for time in range(10,10+1,1):

	# plot name
	csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])
	plotname=''.join([plot_loc,'/',str(time),'.png'])
	plotname_pdf=''.join([plot_loc,'/',str(time),'.pdf'])

	# get dimensional time
	stats_line_num = num_header_lines + (time * model_output_dt) 
	f=open(stats_file)
	line=f.readlines()[stats_line_num]
	time_dim=float(line.split()[1])/1.e6 # Ma

	print("%.0f: t = %.1f Ma" % (time,time_dim))
	model_data  = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

	# print "interpolating model outputs to regular grid..."
	visc   = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,visc_col], (X_low, Y_low), method='linear')
	llith  = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,c_llith_col], (X_low2, Y_low2), method='linear')
	ulith  = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,c_ulith_col], (X_low2, Y_low2), method='linear')

	# get llith contour 
	comp_contour_val = 0.5
	llith_cont = plt.contour(X_low2/1.e3, (ymax-Y_low2)/1.e3, llith, levels=[comp_contour_val])
	ulith_cont = plt.contour(X_low2/1.e3, (ymax-Y_low2)/1.e3, ulith, levels=[comp_contour_val])
	llith_points_tmp = llith_cont.collections[0].get_paths()[0].vertices
	ulith_points_tmp = ulith_cont.collections[0].get_paths()[0].vertices
	# remove shallow and deep points
	cutoff_shall = 110.; cutoff_deep  = 575.; n = 0
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
	# find along-slab distance corresponding to analysis depth
	dz = 1.e9
	for i in range(1,len(llith_points)):
		dz_tmp = np.abs(llith_points[i,1]-(analysis_depth/1.e3))
		if dz_tmp < dz:
			dz = dz_tmp
			analysis_ind = i
	analysis_depth_new = llith_points[analysis_ind,1]
	analysis_slab_dist = llith_points[analysis_ind,2]
	print("pulling out shear stress @ z = %.1f and ds = %.1f km" % (analysis_depth_new,analysis_slab_dist))

	# extract stresses along mid-slab 
	sxx_slab = griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,sxx_col], (llith_points[:,0], llith_points[:,1]), method='linear')
	sxy_slab = griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,sxy_col], (llith_points[:,0], llith_points[:,1]), method='linear')
	# vx_slab =  griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,vx_col],  (llith_points[:,0], llith_points[:,1]), method='linear')
	# vy_slab =  griddata((model_data[:,x_col]/1.e3, (ymax-model_data[:,y_col])/1.e3), model_data[:,vy_col],  (llith_points[:,0], llith_points[:,1]), method='linear')

	# get dip of surface and compute slab-norm velocity
	dips = np.zeros((len(llith_points),1))
	for i in range(len(llith_points)):
		if i == 0:
			dx = llith_points[i+1,0] - llith_points[i,0]
			dy = llith_points[i+1,1] - llith_points[i,1]   
		elif i == len(llith_points):
			dx = llith_points[i+1,0] - llith_points[i-1,0]
			dy = llith_points[i+1,1] - llith_points[i-1,1]
		else:
			dx = llith_points[i,0]   - llith_points[i-1,0]
			dy = llith_points[i,1]   - llith_points[i-1,1] 

		if dx < 0:
			dx = -1. * dx
			dips[i] = 180. - np.rad2deg(np.arctan(dy/dx))
		else:
			dips[i] = np.rad2deg(np.arctan(dy/dx))
	dips=savgol_filter(dips[:,0],251,3)

	# get shear stress perpendicular to slab
	slab_shear_stress = np.zeros((len(llith_points),1))
	for i in range(0,len(sxx_slab)):
		norm_contrib = 2.*sxx_slab[i]*np.sin(np.deg2rad(dips[i]))*np.cos(np.deg2rad(dips[i]));
		shear_contrib = sxy_slab[i]*((np.cos(np.deg2rad(dips[i])))**2 - (np.sin(np.deg2rad(dips[i])))**2);
		slab_shear_stress[i] = norm_contrib - shear_contrib
	slab_shear_stress_splinefit  = splrep(llith_points[:,2],slab_shear_stress,k=3,s=0)
	slab_shear_stress_splinevals = splev(llith_points[:,2],slab_shear_stress_splinefit)
	slab_shear_stress_splinevals = savgol_filter(slab_shear_stress_splinevals,101,3)

	# get down-slab gradient in shear stress
	gradient_slab_shear_stress = np.gradient(slab_shear_stress_splinevals, llith_points[:,2]*1.e3)

	##### get DP and dip #####################
	analysis_depth_km = analysis_depth/1.e3
	dz = 2.0e3
	# get horiz profile in mid-mantle 
	midmant_prof_loc = ymax - analysis_depth                          
	midmant_prof = model_data[model_data[:,y_col] < (midmant_prof_loc+dz)] 
	midmant_prof = midmant_prof[midmant_prof[:,y_col] > (midmant_prof_loc-dz)] 
	midmant_prof = midmant_prof[midmant_prof[:,x_col].argsort()] # sort by x
	# get slab location
	# left (sub-slab):
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
	print("x_left = %.2f (%.1f km depth), x_right = %.2f km (%.1f km depth)," % (x_left/1.e3,y_left,x_right/1.e3,y_right))
	print("x_center = %.2f (%.1f km depth)" % (x_center/1.e3,y_center))
	# dip
	dz = 1e9;
	for i in range(1,len(llith_points)):
	        dz_tmp = np.abs(llith_points[i,1]-analysis_depth_km)
	        if dz_tmp < dz:
	                dz = dz_tmp
	                dip_ind = i
	                dip_midmant=dips[i]
	print("midmant dip = %.2f deg" % dip_midmant)
	# slab  normal thickness
	slabnorm_thick = (x_right-x_left) * np.sin(np.deg2rad(dip_midmant)) # m
	print("midmant slab thick  = %.2f km" % (slabnorm_thick/1.e3))
	# compute theoretical DP
	DP_anal = drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(dip_midmant)) # Pa


	# pull out DP
	z_shallow = y_center - (0.5*slabnorm_thick*1.e-3)*np.cos(np.deg2rad(dip_midmant)) # km
	z_deep    = y_center + (0.5*slabnorm_thick*1.e-3)*np.cos(np.deg2rad(dip_midmant)) # km
	shallow_prof_loc = ymax - (z_shallow*1.e3)                     
	dz = 2.0e3
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
	    y = (ymax - shallow_prof[i,y_col])/1.e3
	    c_ulith = shallow_prof[i,c_ulith_col]
	    if x > slab_x_right and c_ulith > 0.5:
	        slab_x_right = x
	        slab_y_right = y
	        ind_right = i
	slab_x_left = 1e9;
	for i in range(len(deep_prof)):
	    x = deep_prof[i,x_col]
	    y = (ymax - deep_prof[i,y_col])/1.e3
	    c_llith = deep_prof[i,c_llith_col]
	    if x < slab_x_left and c_llith > 0.5:
	        slab_x_left = x
	        slab_y_left = y
	        ind_left = i
	# find pressure at these points
	ds = 10.e3
	misfit_left = 1.e9
	misfit_right = 1.e9
	pos_left = slab_x_left - ds
	pos_right = slab_x_right + ds
	for i in range(len(shallow_prof)):
	    x = shallow_prof[i,x_col]
	    y = (ymax - shallow_prof[i,y_col])/1.e3
	    misfit_right_tmp = np.abs(x-pos_right)
	    if misfit_right_tmp < misfit_right:
	    	misfit_right = misfit_right_tmp
	    	Pright       = shallow_prof[i,P_col] 
	    	Pright_x     = shallow_prof[i,x_col]
	for i in range(len(deep_prof)):
	    x = deep_prof[i,x_col]
	    y = (ymax - deep_prof[i,y_col])/1.e3
	    misfit_left_tmp  = np.abs(x-pos_left)
	    if misfit_left_tmp < misfit_left:
	    	misfit_left = misfit_left_tmp
	    	Pleft       = deep_prof[i,P_col] 
	    	Pleft_x     = deep_prof[i,x_col]
	DP_mod = Pleft - Pright # Pa
	DP_missing = DP_anal - DP_mod

	################# plotting visc, density fields #######################
	fig=plt.figure()
	gs=GridSpec(5,1) 

	# plot viscosity field
	ax1=fig.add_subplot(gs[0,0])
	visc_plot = ax1.contourf(X_low/1.e3, (ymax-Y_low)/1.e3, np.log10(visc), cmap=cm.get_cmap('plasma_r'),levels=np.linspace(19,24,501))
	ax1.set_ylim([(ymax-ymin_plot)/1.e3,0])   
	ax1.tick_params(direction='out',length=2, labelsize=6)
	ax1.annotate(''.join([str("%.1f" % (time_dim)),' Myr']), xy=(0.01,0.6), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=10,color='k')		
	ax1.annotate(''.join(['DP: mod = ',str("%.1f" % (DP_mod/1.e6)),', anal = ',str("%.1f" % (DP_anal/1.e6)),', missing =',str("%.1f" % (DP_missing/1.e6)),' MPa']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')		
	ax1.plot(llith_points[:,0], llith_points[:,1], linewidth=0.5, color='black',zorder=6)
	# ax1.plot([((x_left)/1.e3)-125, ((x_right)/1.e3)+125.],[analysis_depth_km,analysis_depth_km],linewidth=0.5, color='red',zorder=6)
	ax1.scatter(slab_x_right/1.e3, slab_y_right, s=0.5,color='red',zorder=3)
	ax1.scatter(slab_x_left/1.e3, slab_y_left, s=0.5,color='red',zorder=3)
	ax1.axis('equal')

	# # plot DP extraction gradient
	# ax2=fig.add_subplot(gs[1,0])
	# ax2.set_ylabel("P [MPa]",size=6.5)
	# ax2.set_xlim(((slab_x_left)/1.e3)-125.,((slab_x_left)/1.e3)+125.); 
	# ax2.set_ylim(-50,50); 
	# ax2.tick_params(axis='x', labelsize=6)
	# ax2.tick_params(axis='y', labelsize=6)
	# ax2.plot(shallow_prof[:,x_col]/1.e3,shallow_prof[:,P_col]/1.e6,color='blue', linestyle='-', linewidth=0.75, zorder=2) 
	# ax2.plot(deep_prof[:,x_col]/1.e3,deep_prof[:,P_col]/1.e6,color='black', linestyle='-', linewidth=0.75, zorder=2) 
	# ax2.scatter(Pleft_x/1.e3, Pleft/1.e6, s=0.5,color='black',zorder=3)
	# ax2.scatter(Pright_x/1.e3,Pright/1.e6,s=0.5,color='blue',zorder=3)
	# ax2.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)
	# ax2.axvline(x=slab_x_left/1.e3, color='red',linestyle='-',linewidth=0.75, zorder=0)
	# ax2.axvline(x=slab_x_right/1.e3, color='red',linestyle='-',linewidth=0.75, zorder=0)
	# ax2.annotate(''.join(['DP: mod = ',str("%.1f" % (DP_mod/1.e6)),', anal = ',str("%.1f" % (DP_anal/1.e6)),', missing =',str("%.1f" % (DP_missing/1.e6)),' MPa']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')		

	# plot mid-slab stresses
	ax3=fig.add_subplot(gs[1,0])
	ax3.set_ylabel(r"$\sigma$ [MPa]",size=6.5)
	ax3.set_xlim(0,600); 
	ax3.set_ylim(-60,60);
	ax3.tick_params(axis='x', labelsize=6)
	ax3.tick_params(axis='y', labelsize=6)
	ax3.plot(llith_points[:,2],sxx_slab/1.e6,label='sxx',color='darkblue', linestyle='-', linewidth=1.5, zorder=2)   
	ax3.plot(llith_points[:,2],sxy_slab/1.e6,label='sxy',color='red',      linestyle='-', linewidth=1.5, zorder=2)
	ax3.plot(llith_points[:,2],-1.*sxx_slab/1.e6,label='syy',color='forestgreen', linestyle='-', linewidth=1.5, zorder=2)
	ax3.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)
	ax3.axvline(x=analysis_slab_dist, color='peachpuff',linestyle='-',linewidth=4, zorder=0)
	ax3.legend(fontsize=6,loc='right',facecolor='white',fancybox=True, framealpha=1)

	# plot dip
	ax4=fig.add_subplot(gs[2,0])
	ax4.set_ylabel(r"$\theta$ [deg]",size=6.5)
	ax4.set_xlim(0,600); 
	ax4.set_ylim(30,100); 
	ax4.tick_params(axis='x', labelsize=6)
	ax4.tick_params(axis='y', labelsize=6)
	ax4.plot(llith_points[:,2],dips,color='black', linestyle='-', linewidth=1.5)  
	ax4.axhline(y=90, color='gray',linestyle='--',linewidth=1, zorder=1)
	ax4.axvline(x=analysis_slab_dist, color='peachpuff',linestyle='-',linewidth=4, zorder=0)

	# plot stress
	ax5=fig.add_subplot(gs[3,0])
	ax5.set_ylabel(r"$\tau_{n}$ [MPa]",size=6.5)
	ax5.set_xlim(0,600); 
	ax5.set_ylim(0,70); 
	ax5.tick_params(axis='x', labelsize=6)
	ax5.tick_params(axis='y', labelsize=6)
	ax5.plot(llith_points[:,2],slab_shear_stress_splinevals/1.e6,color='black', linestyle='-', linewidth=1.5,zorder=2)  
	ax5.axvline(x=analysis_slab_dist, color='peachpuff',linestyle='-',linewidth=4, zorder=0)

	# plot stress gradient
	ax6=fig.add_subplot(gs[4,0])
	ax6.set_ylabel(r"h*d$\tau_{n}$/ds   [MPa]",size=6.5)
	ax6.set_xlabel('along-slab distance [km]',size=6.5)
	ax6.set_xlim(0,600); 
	ax6.set_ylim(-20,20); 
	ax6.tick_params(axis='x', labelsize=6)
	ax6.tick_params(axis='y', labelsize=6)
	ax6.plot(llith_points[:,2],slabnorm_thick*(gradient_slab_shear_stress/1.e6),color='black', linestyle='-', linewidth=1.5, zorder=2)   
	ax6.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)
	ax6.axvline(x=analysis_slab_dist, color='peachpuff',linestyle='-',linewidth=4, zorder=0)

	print("saving fields figure to %s..." % plotname)
	plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
	# plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
	plt.clf()

	 