#!/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
from functions import create_grid, get_slab_midplane, extract_horiz_prof, get_slablocation_from_horiz_prof

model_name=str(sys.argv[1])            
time=int(sys.argv[2])  

analysis_depth1 = 230e3
analysis_depth2 = 330e3
analysis_depth3 = 430e3

# model properties
xmax=5800.e3
ymax=1450.e3

# ASPECT output 
csvs_loc =  'csv_outputs/'
models_loc =  'raw_outputs/'
stats_file = ''.join([models_loc,str(model_name),'/statistics'])
model_output_dt  = 50 # output dt as set in ASPECT .prm file (for getting the dimensional time)
num_header_lines = 16 # num header lines in stats_files (for getting the dimensional time)

# where to put the plots
plot_loc   = ''.join(['plots/evolution/simple/zoomed/',str(model_name)])
if not os.path.exists(plot_loc):
    os.mkdir(plot_loc)

# column numbers of the relevant properties in the .csv file. 
c_crust_col=23; c_ulith_col = 24;   c_llith_col = 25;
P_col = 29;     x_col=30;           y_col = 31;
visc_col=26;    vx_col=0;           vz_col=1;

# create low res grid for plotting
xmin_plot=2500e3; xmax_plot = 4200e3
ymin_plot=0; grid_res=2.0e3
X_low, Y_low = create_grid(xmin_plot,xmax_plot,ymin_plot,ymax,grid_res)
# extremely low res grid for plotting velocities
X_vels, Y_vels = create_grid(xmin_plot,xmax_plot,ymin_plot,ymax,10.e3)
# higher res grid for detailed calculations
xmin_plot2=2500.e3; xmax_plot2=3750.e3
ymin_plot2=ymax-600.e3; grid_res2=1.0e3
X_low2, Y_low2 = create_grid(xmin_plot,xmax_plot,ymin_plot2,ymax,grid_res2)

# plot name
csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])
plotname=''.join([plot_loc,'/',str(time),'.png'])
plotname_pdf=''.join([plot_loc,'/',str(time),'.pdf'])

# get dimensional time
stats_line_num = num_header_lines + (time * model_output_dt) 
f=open(stats_file); line=f.readlines()[stats_line_num]
time_dim=float(line.split()[1])/1.e6 # Myr
print("%.0f: t = %.1f Myr" % (time,time_dim))
model_data  = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

# print "interpolating model outputs to regular grid..."
visc   = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,visc_col], (X_low, Y_low), method='nearest')
llith  = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,c_llith_col], (X_low2, Y_low2), method='linear')
vx     = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,vx_col],   (X_vels, Y_vels), method='nearest')
vz     = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,vz_col],   (X_vels, Y_vels), method='nearest')
vmag   = np.sqrt(vx**2 + vz**2)

# get lithosphere contour, and then trim it to get slab mid-plane
comp_contour_val = 0.5
llith_cont = plt.contour(X_low2/1.e3, (ymax-Y_low2)/1.e3, llith, levels=[comp_contour_val])
llith_points_tmp = llith_cont.collections[0].get_paths()[0].vertices
cutoff_shall = 110.; cutoff_deep  = 575.;
llith_points = get_slab_midplane(llith_points_tmp,cutoff_shall,cutoff_deep)

# get mid-point at the analysis depth
misfit = 1e9
for j in range(len(llith_points)):
    misfit_tmp = np.abs(llith_points[j,1]-(analysis_depth1/1.e3))
    if misfit_tmp < misfit:
        misfit = misfit_tmp
        xcent1 = llith_points[j,0]
        zcent1 = llith_points[j,1]
# get mid-point at the analysis depth
misfit = 1e9
for j in range(len(llith_points)):
    misfit_tmp = np.abs(llith_points[j,1]-(analysis_depth2/1.e3))
    if misfit_tmp < misfit:
        misfit = misfit_tmp
        xcent2 = llith_points[j,0]
        zcent2 = llith_points[j,1]
# get mid-point at the analysis depth
misfit = 1e9
for j in range(len(llith_points)):
    misfit_tmp = np.abs(llith_points[j,1]-(analysis_depth3/1.e3))
    if misfit_tmp < misfit:
        misfit = misfit_tmp
        xcent3 = llith_points[j,0]
        zcent3 = llith_points[j,1]


###################### plotting #########################
fig=plt.figure()
gs=GridSpec(3,1) 

# plot viscosity field
ax1=fig.add_subplot(gs[0,0])
visc_plot = ax1.contourf(X_low/1.e3, (ymax-Y_low)/1.e3, np.log10(visc), cmap=cm.get_cmap('hot_r'),levels=np.linspace(20,25,501),extend='max')
# flow vectors
ax1.streamplot(np.flipud(X_vels/1.e3),np.flipud((ymax-Y_vels)/1.e3),np.flipud(vx),np.flipud(-1.0*vz),color='dimgray', linewidth=vmag*18, arrowsize=0.3,density=1.2)
# color bar
cbar = plt.colorbar(visc_plot, cax = fig.add_axes([0.92, 0.7, 0.015, 0.15]), ticks=[20,21,22,23,24,25], ticklocation = 'right')
cbar.ax.tick_params(axis='y',labelsize=5.5,pad=1,left=False,labelleft=False,right=True,labelright=True)
# axis properties
ax1.axis('equal')
ax1.set_ylim([1150,0])   
ax1.set_xlim([2500,4200])
ax1.tick_params(direction='out',length=2, labelsize=6)
ax1.annotate(''.join([str("%.1f" % (time_dim)),' Myr']), xy=(0.01,0.5), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=10,color='k') 
ax1.scatter(xcent1,  zcent1,  s=3,color='gray', zorder=3)
ax1.scatter(xcent2,  zcent2,  s=3,color='black',zorder=3)
ax1.scatter(xcent3,  zcent3,  s=3,color='gray', zorder=3)


# ax1.plot(llith_points[:,0], llith_points[:,1], linewidth=0.5, color='black',zorder=6)

# print("saving fields figure to %s..." % plotname)
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()
