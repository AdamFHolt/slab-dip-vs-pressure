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
from scipy.interpolate import splrep, splev

model_name=str(sys.argv[1])	            
time=int(sys.argv[2])   # largest number in csv_outputs/ filenames

xmax=5800.e3
ymax=1450.e3

mant_visc = 2.5e20
slab_visc = mant_visc * 500.

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
# Note: these numberings will change as outputted model properties (in .prm file) change.
visc_col=26; 	vx_col=0;  		 	vy_col=1;
c_crust_col=23; c_ulith_col = 24; 	c_llith_col = 25;
P_col = 29; 	x_col=30;  		 	y_col = 31;
sxx_col = 3;	syy_col = 7;		sxy_col = 4; 

# plot name
csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])
plotname=''.join([plot_loc,'/',str(time),'.png'])
plotname_pdf=''.join([plot_loc,'/',str(time),'.pdf'])

# get dimensional time
stats_line_num = num_header_lines + (time * model_output_dt) 
f=open(stats_file)
line=f.readlines()[stats_line_num]
time_dim=float(line.split()[1])/1.e6 # Ma

print "%.0f: t = %.1f Ma" % (time,time_dim)
model_data  = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

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

# print "interpolating model outputs to regular grid..."
visc   = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,visc_col], (X_low, Y_low), method='linear')
llith  = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,c_llith_col], (X_low2, Y_low2), method='linear')
# sxx    = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,sxx_col], (X_low2, Y_low2), method='linear')
# sxy    = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,sxy_col], (X_low2, Y_low2), method='linear')

# get llith contour 
comp_contour_val = 0.5
llith_cont = plt.contour(X_low2/1.e3, (ymax-Y_low2)/1.e3, llith, levels=[comp_contour_val],linewidths=0.5, colors='red',zorder=4)   
llith_points_tmp = llith_cont.collections[0].get_paths()[0].vertices
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
slab_shear_stress_splinefit  = splrep(llith_points[:,2],slab_shear_stress[:,0],k=5,s=3)
slab_shear_stress_splinevals = splev(llith_points[:,2],slab_shear_stress_splinefit)/1.e6

# get down-slab gradient in shear stress
gradient_slab_shear_stress = np.gradient(slab_shear_stress[:,0], llith_points[:,2]*1.e3)


##### get DP and dip #####################
midmant_depth = 330.e3
# get horiz profile in mid-mantle 
midmant_prof_loc = ymax - midmant_depth                          
midmant_prof = model_data[model_data[:,y_col] < (midmant_prof_loc+3.0e3)] 
midmant_prof = midmant_prof[midmant_prof[:,y_col] > (midmant_prof_loc-3.0e3)] 
midmant_prof = midmant_prof[midmant_prof[:,x_col].argsort()] # sort by x
# compute dip from llith contour
# x_shall = 0; x_deep = 0;
# for d in range(len(tc.collections[0].get_paths())): 
# 	p = tc.collections[0].get_paths()[d]
# 	x = p.vertices[:,0]; z = p.vertices[:,1]
# 	for j in range(len(x)):
# 		if x[j] > x_shall and z[j] < (235.+(grid_res/1.e3)) and z[j] > (235.-(grid_res/1.e3)):
# 			x_shall = x[j]; z_shall = z[j]
# 		if x[j] > x_deep and  z[j] < (265.+(grid_res/1.e3)) and z[j] > (265.-(grid_res/1.e3)):
# 			x_deep = x[j];  z_deep = z[j]
#     if x_deep > 0:
#         dip_deep = np.rad2deg(np.arctan((z_deep-z_shall)/(x_deep-x_shall)))
# else:
# 	dip_deep = 0


################# plotting visc, density fields #######################
fig=plt.figure()
gs=GridSpec(5,1) 

# plot viscosity field
ax1=fig.add_subplot(gs[0,0])
visc_plot = ax1.contourf(X_low/1.e3, (ymax-Y_low)/1.e3, np.log10(visc), cmap=cm.get_cmap('plasma_r'),levels=np.linspace(19,24,501))
ax1.set_ylim([(ymax-ymin_plot)/1.e3,0])   
ax1.tick_params(direction='out',length=2, labelsize=6)
ax1.annotate(''.join([str("%.1f" % (time_dim)),' Myr']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=10,color='k')		
ax1.plot(llith_points[:,0], llith_points[:,1], linewidth=0.5, color='black',zorder=6)

# plot mid-slab stresses
ax3=fig.add_subplot(gs[1,0])
ax3.set_ylabel('stress [MPa]',size=6.5)
ax3.set_xlim(0,600); 
ax3.set_ylim(-60,60);
ax3.tick_params(axis='x', labelsize=6)
ax3.tick_params(axis='y', labelsize=6)
ax3.plot(llith_points[:,2],sxx_slab/1.e6,label='sxx',color='darkblue', linestyle='-', linewidth=1.5)   
ax3.plot(llith_points[:,2],sxy_slab/1.e6,label='sxy',color='red',      linestyle='-', linewidth=1.5)
ax3.plot(llith_points[:,2],-1.*sxx_slab/1.e6,label='syy',color='forestgreen', linestyle='-', linewidth=1.5)   
ax3.legend(fontsize=6,loc='right',facecolor='white',fancybox=True, framealpha=1)

# plot dip
ax4=fig.add_subplot(gs[2,0])
ax4.set_ylabel('dip [degrees]',size=6.5)
# ax4.set_xlabel('along-slab distance',size=6.5)
ax4.set_xlim(0,600); 
ax4.set_ylim(30,100); 
ax4.tick_params(axis='x', labelsize=6)
ax4.tick_params(axis='y', labelsize=6)
ax4.plot(llith_points[:,2],dips,color='black', linestyle='-', linewidth=1.5)  

# plot stress
ax5=fig.add_subplot(gs[3,0])
ax5.set_ylabel('shear stress [MPa]',size=6.5)
ax5.set_xlim(0,600); 
ax5.set_ylim(0,70); 
ax5.tick_params(axis='x', labelsize=6)
ax5.tick_params(axis='y', labelsize=6)
ax5.plot(llith_points[:,2],slab_shear_stress_splinevals,color='black', linestyle='-', linewidth=1.5,zorder=2)  

# plot stress gradient
ax5=fig.add_subplot(gs[4,0])
ax5.set_ylabel('h*(dtau/ds) [MPa]',size=6.5)
ax5.set_xlabel('along-slab distance [km]',size=6.5)
ax5.set_xlim(0,600); 
ax5.set_ylim(-100,100); 
ax5.tick_params(axis='x', labelsize=6)
ax5.tick_params(axis='y', labelsize=6)
ax5.plot(llith_points[:,2],100.e3*(gradient_slab_shear_stress/1.e6),color='black', linestyle='-', linewidth=1.5)   

print("saving fields figure to %s..." % plotname)
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()

  
