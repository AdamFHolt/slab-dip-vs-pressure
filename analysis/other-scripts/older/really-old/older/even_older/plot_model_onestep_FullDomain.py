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
from functions import create_grid, get_slab_midplane, get_dip_slab_midplane
from functions import extract_horiz_prof, get_slablocation_from_horiz_prof
from functions import get_nearslab_pressures, get_dip_at_certain_depth
from functions import get_stress_profile, convert_to_slabnorm_shearstress

model_name=str(sys.argv[1])            
max_time=int(sys.argv[2])   # largest number in csv_outputs/ filenames
analysis_depth = float(sys.argv[3])     # m (depth for DP extraction and central point of shear stress derivative)
analysis_depth_dz = float(sys.argv[4])  # m (depth interval for shear stress derivative)
ds = float(sys.argv[5])                 # m (distance from slab to pull out DP)
dz = float(sys.argv[6])                 # m (height used to extract horizontal profiles, i.e., points +/- this dz)


# analysis_depth = 330.e3   # m (depth for DP extraction and central point of shear stress derivative)
# analysis_depth_dz = 15.e3     # m (depth interval for shear stress derivative)
# ds = 7.5e3                    # m (distance from slab to pull out DP)
# dz = 1.0e3                    # m (height used to extract horizontal profiles, i.e., points +/- this dz)

# model properties
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
saved_stresses_name = ''.join(['text_files/lo-res/',model_name,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])

# where to put the plots
plot_loc   = ''.join(['plots/evolution/',str(model_name)])
if not os.path.exists(plot_loc):
    os.mkdir(plot_loc)

# column numbers of the relevant properties in the .csv file. 
visc_col=26;    vx_col=0;           vy_col=1;
c_crust_col=23; c_ulith_col = 24;   c_llith_col = 25;
P_col = 29;     x_col=30;           y_col = 31;
sxx_col = 3;    syy_col = 7;        sxy_col = 4; 

# create low res grid for plotting
xmin_plot=0; ymin_plot=500.e3; grid_res=2.5e3
X_low, Y_low = create_grid(xmin_plot,xmax,ymin_plot,ymax,grid_res)

# higher res grid for detailed calculations
xmin_plot2=2500.e3; xmax_plot2=3750.e3
ymin_plot2=ymax-600.e3; grid_res2=0.25e3
X_low2, Y_low2 = create_grid(xmin_plot2,xmax_plot2,ymin_plot2,ymax,grid_res2)

first_time=8
saved_stresses = np.zeros(((max_time-first_time)+1,10)) # time [Myr], DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip, stress term, stress term B, stress term C,slabnorm thick
ind = 0 

for time in range(first_time,max_time+1,1):
# for time in np.array([15]):

    # plot name
    csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])
    plotname=''.join([plot_loc,'/',str(time),'.png'])
    plotname_pdf=''.join([plot_loc,'/',str(time),'.pdf'])

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
    dips = savgol_filter(dips[:,0],151,3)

    # get dips at the relevant depths
    dip_midmant       = get_dip_at_certain_depth(dips,llith_points,analysis_depth)
    dip_midmant_shall = get_dip_at_certain_depth(dips,llith_points,analysis_depth_shall)
    dip_midmant_deep  = get_dip_at_certain_depth(dips,llith_points,analysis_depth_deep)
    # print("dips: shallow = %.2f, center = %.2f, deep = %.2f deg" % (dip_midmant_shall,dip_midmant,dip_midmant_deep))

    # slab normal thicknesses at the relevant depths
    slabnorm_thick       = (x_right-x_left) * np.sin(np.deg2rad(dip_midmant))                   # m
    slabnorm_thick_shall = (x_right_shall-x_left_shall) * np.sin(np.deg2rad(dip_midmant_shall)) # m
    slabnorm_thick_deep  = (x_right_deep -x_left_deep)  * np.sin(np.deg2rad(dip_midmant_deep))  # m

    # compute theoretical DP
    DP_anal = drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(dip_midmant)) # Pa

    # pull out model DP
    Pleft, Pright, Pleft_x, Pright_x, Pleft_y, Pright_y, slab_x_left, slab_y_left, slab_x_right, slab_y_right = \
        get_nearslab_pressures(y_center,slabnorm_thick,dip_midmant,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,dz) # Pa
    DP_mod = Pleft - Pright # Pa
    DP_missing = DP_anal - DP_mod
    print("DPs: analytical = %.1f, model %.1f, missing = %.1f MPa" % (DP_anal/1.e6,DP_mod/1.e6,DP_missing/1.e6))

    # other depths for shear stress
    Pleft_shall, Pright_shall, Pleft_x_shall, Pright_x_shall, Pleft_y_shall, Pright_y_shall, slab_x_left_shall, slab_y_left_shall, slab_x_right_shall, slab_y_right_shall = \
        get_nearslab_pressures(y_center_shall,slabnorm_thick_shall,dip_midmant_shall,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,dz) # Pa
    Pleft_deep, Pright_deep, Pleft_x_deep, Pright_x_deep, Pleft_y_deep, Pright_y_deep, slab_x_left_deep, slab_y_left_deep, slab_x_right_deep, slab_y_right_deep = \
        get_nearslab_pressures(y_center_deep, slabnorm_thick_deep, dip_midmant_deep, model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,dz) # Pa  
    DP_mod_shall = Pleft_shall - Pright_shall
    DP_mod_deep  = Pleft_deep  - Pright_deep

    # get slab stresses along slab-perp profiles
    sxx_slab, sxy_slab, profile = \
        get_stress_profile(model_data,slab_x_left,slab_y_left,slab_x_right,slab_y_right,ymax,x_col,y_col,sxx_col,sxy_col)
    sxx_slab_shall, sxy_slab_shall, profile_shall = \
        get_stress_profile(model_data,slab_x_left_shall,slab_y_left_shall,slab_x_right_shall,slab_y_right_shall,ymax,x_col,y_col,sxx_col,sxy_col)
    sxx_slab_deep, sxy_slab_deep, profile_deep = \
        get_stress_profile(model_data,slab_x_left_deep,slab_y_left_deep,slab_x_right_deep,slab_y_right_deep,ymax,x_col,y_col,sxx_col,sxy_col)

    # rotate to get shear stresses
    slab_shear_stress       = convert_to_slabnorm_shearstress(sxx_slab,      sxy_slab,      dip_midmant)
    slab_shear_stress_shall = convert_to_slabnorm_shearstress(sxx_slab_shall,sxy_slab_shall,dip_midmant_shall)
    slab_shear_stress_deep  = convert_to_slabnorm_shearstress(sxx_slab_deep, sxy_slab_deep, dip_midmant_deep)

    # cut out the boundaries of the slabs (which have stress discontinuities)
    inds                        = np.where((profile[:,2]>7) & (profile[:,2]<profile[len(profile)-1,2]-5.))[0]
    slab_shear_stress_cut       = slab_shear_stress[inds]
    profile_cut                 = profile[inds,:]
    inds_shall                  = np.where((profile_shall[:,2]>7) & (profile_shall[:,2]<profile_shall[len(profile_shall)-1,2]-5.))[0]
    slab_shear_stress_shall_cut = slab_shear_stress_shall[inds_shall]
    profile_shall_cut           = profile_shall[inds_shall,:]
    inds_deep                   = np.where((profile_deep[:,2]>7) & (profile_deep[:,2]<profile_deep[len(profile_deep)-1,2]-5.))[0]
    slab_shear_stress_deep_cut  = slab_shear_stress_deep[inds_deep]
    profile_deep_cut            = profile_deep[inds_deep,:]

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
    print("shear stresses: shallow = %.1f; center = %.1f; deep = %.1f MPa" % (slab_stress_shall/1.e6,slab_stress/1.e6,slab_stress_deep/1.e6))

    # get d(stress)/ds using central difference approximation
    ds_shall = np.sqrt((x_center - x_center_shall)**2 + (y_center*1e3 - y_center_shall*1e3)**2)
    ds_deep  = np.sqrt((x_center - x_center_deep)**2  + (y_center*1e3 - y_center_deep*1e3)**2)

    # # try getting more accurate distances from mid-plane
    # misfit_shall  = 1.e9
    # misfit_center = 1.e9
    # misfit_deep   = 1.e9
    # for i in range(len(llith_points)):
    # 	misfit_shall_tmp  = np.sqrt((x_center_shall - llith_points[i,0]*1e3)**2  + (y_center_shall*1e3 - llith_points[i,1]*1e3)**2)
    # 	misfit_center_tmp = np.sqrt((x_center       - llith_points[i,0]*1e3)**2  + (y_center*1e3       - llith_points[i,1]*1e3)**2)
    # 	misfit_deep_tmp	  = np.sqrt((x_center_deep  - llith_points[i,0]*1e3)**2  + (y_center_deep*1e3  - llith_points[i,1]*1e3)**2)
    # 	if misfit_shall_tmp < misfit_shall:
    # 		misfit_shall  = misfit_shall_tmp
    # 		shall_ind = i
    # 	if misfit_center_tmp < misfit_center:
    # 		misfit_center  = misfit_center_tmp
    # 		center_ind = i
    # 	if misfit_deep_tmp < misfit_deep:
    # 		misfit_deep  = misfit_deep_tmp
    # 		deep_ind = i

    # print(shall_ind,center_ind,deep_ind)
    # print("recalculating distances:")
    # print("shallow: x=%.1f, y=%.1f; x=%.1f, y=%.1f" % (x_center_shall,y_center_shall,llith_points[shall_ind,0],llith_points[shall_ind,1]))
    # print("center: x=%.1f, y=%.1f; x=%.1f, y=%.1f" % (x_center,y_center,llith_points[center_ind,0],llith_points[center_ind,1]))
    # print("deep: x=%.1f, y=%.1f; x=%.1f, y=%.1f" % (x_center_deep,y_center_deep,llith_points[deep_ind,0],llith_points[deep_ind,1]))

    # ds_shall_new = 0
    # for j in range(shall_ind,center_ind,1):
    # 	dist_tmp = np.sqrt((llith_points[j,0]-llith_points[j+1,0])**2  + (llith_points[j,1]-llith_points[j+1,1])**2)
    # 	ds_shall_new = ds_shall_new + dist_tmp
    # ds_deep_new = 0
    # for j in range(center_ind,deep_ind+1,1):
    # 	dist_tmp = np.sqrt((llith_points[j,0]-llith_points[j+1,0])**2  + (llith_points[j,1]-llith_points[j+1,1])**2)
    # 	ds_deep_new = ds_deep_new + dist_tmp
    # print(ds_shall_new,ds_deep_new)

    slab_stress_gradient = (slab_stress_deep - slab_stress_shall)/(ds_shall + ds_deep)
    slab_stress_term = slab_stress_gradient * 0.5*(slabnorm_thick_shall + slabnorm_thick_deep)
    # print(ds_shall/1.e3,ds_deep/1.e3,(ds_shall + ds_deep)/1.e3)
    # print(slabnorm_thick_shall/1.e3,slabnorm_thick_shall/1.e3)
    # alternative ways of calculating the derivatives
    slab_stress_gradient_b = (slab_stress_deep - slab_stress)/(ds_deep)
    slab_stress_term_b = slab_stress_gradient_b * 0.5*(slabnorm_thick + slabnorm_thick_deep)
    slab_stress_gradient_c = (slab_stress - slab_stress_shall)/(ds_shall)
    slab_stress_term_c = slab_stress_gradient_c * 0.5*(slabnorm_thick + slabnorm_thick_shall)
    print("shear stress term = %.2f MPa (others: %.2f, %.2f)" % (slab_stress_term/1.e6,slab_stress_term_b/1.e6,slab_stress_term_c/1.e6))

    saved_stresses[ind,:] = time_dim, DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip_midmant, slab_stress_term, slab_stress_term_b, slab_stress_term_c, slabnorm_thick
    ind = ind + 1

    ###################### plotting #########################
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
    ax1.plot(profile_shall[:,0],profile_shall[:,1], linewidth=0.5, color='green',zorder=6)
    ax1.plot(profile_deep[:,0], profile_deep[:,1],  linewidth=0.5, color='blue',zorder=6)
    ax1.plot(profile[:,0],      profile[:,1],       linewidth=0.5, color='red',zorder=6)
    ax1.scatter(Pleft_x/1.e3,  Pleft_y/1.e3,  s=0.25,color='black',zorder=3)
    ax1.scatter(Pright_x/1.e3, Pright_y/1.e3, s=0.25,color='black',zorder=3)
    #ax1.set_ylim([500.,0])   
    #ax1.set_xlim([2500,3750])   
    ax1.axis('equal')

    # plot mid-slab stresses
    ax2=fig.add_subplot(gs[1,0])
    ax2.set_ylabel(r"shall. $\sigma$ [MPa]",size=6.5)
    ax2.set_xlim(-10,100); 
    ax2.tick_params(axis='x', labelsize=6)
    ax2.tick_params(axis='y', labelsize=6)
    ax2.plot(profile[:,2],sxx_slab_shall/1.e6,label='sxx',color='tan', linestyle='-', linewidth=1.5, zorder=2)   
    ax2.plot(profile[:,2],sxy_slab_shall/1.e6,label='sxy',color='indianred',      linestyle='-', linewidth=1.5, zorder=2)
    ax2.plot(profile[:,2],-1.*sxx_slab_shall/1.e6,label='syy',color='plum', linestyle='-', linewidth=1.5, zorder=2)
    ax2.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)
    ax2.legend(fontsize=6,loc='right',facecolor='white',fancybox=True, framealpha=1)

    ax3=fig.add_subplot(gs[2,0])
    ax3.set_ylabel(r"mid. $\sigma$ [MPa]",size=6.5)
    ax3.set_xlim(-10,100);
    ax3.tick_params(axis='x', labelsize=6)
    ax3.tick_params(axis='y', labelsize=6)
    ax3.plot(profile[:,2],sxx_slab/1.e6,label='sxx',color='tan', linestyle='-', linewidth=1.5, zorder=2)
    ax3.plot(profile[:,2],sxy_slab/1.e6,label='sxy',color='indianred',      linestyle='-', linewidth=1.5, zorder=2)
    ax3.plot(profile[:,2],-1.*sxx_slab/1.e6,label='syy',color='plum', linestyle='-', linewidth=1.5, zorder=2)
    ax3.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)

    ax4=fig.add_subplot(gs[3,0])
    ax4.set_ylabel(r"deep $\sigma$ [MPa]",size=6.5)
    ax4.set_xlim(-10,100);
    ax4.tick_params(axis='x', labelsize=6)
    ax4.tick_params(axis='y', labelsize=6)
    ax4.plot(profile[:,2],sxx_slab_deep/1.e6,label='sxx',color='tan', linestyle='-', linewidth=1.5, zorder=2)
    ax4.plot(profile[:,2],sxy_slab_deep/1.e6,label='sxy',color='indianred',      linestyle='-', linewidth=1.5, zorder=2)
    ax4.plot(profile[:,2],-1.*sxx_slab_deep/1.e6,label='syy',color='plum', linestyle='-', linewidth=1.5, zorder=2)
    ax4.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)

    # plot stress
    ax5=fig.add_subplot(gs[4,0])
    ax5.set_ylabel(r"$\tau_{n}$ [MPa]",size=6.5)
    ax5.set_xlim(-10,100); 
    ax5.tick_params(axis='x', labelsize=6)
    ax5.tick_params(axis='y', labelsize=6)
    ax5.plot(profile[:,2],      slab_shear_stress/1.e6,      color='red', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
    ax5.plot(profile_shall[:,2],slab_shear_stress_shall/1.e6, color='green', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
    ax5.plot(profile_deep[:,2], slab_shear_stress_deep/1.e6,  color='blue', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
    ax5.plot(profile_cut[:,2],      slab_shear_stress_cut/1.e6,       color='red', linestyle='-', linewidth=1.5,zorder=3)  
    ax5.plot(profile_shall_cut[:,2],slab_shear_stress_shall_cut/1.e6, color='green', linestyle='-', linewidth=1.5,zorder=3)  
    ax5.plot(profile_deep_cut[:,2], slab_shear_stress_deep_cut/1.e6,  color='blue', linestyle='-', linewidth=1.5,zorder=3) 
    ax5.annotate(''.join(['slab stress term = ',str("%.1f" % (slab_stress_term/1.e6)),' MPa']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=8,color='k')        
    ax5.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)

    # print("saving fields figure to %s..." % plotname)
    plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
    plt.clf()

np.savetxt(saved_stresses_name, saved_stresses)
