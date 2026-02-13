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

def get_vs_near_slab_center(horiz_prof, x_center, y_center_km, ymax, x_col, y_col, vx_col, vy_col, dip_deg):
    """
    Extract slab-parallel velocity at the point in a horizontal profile nearest
    the slab-center location (x_center, y_center_km depth).
    """
    if len(horiz_prof) == 0:
        return np.nan, np.nan, np.nan

    y_center_model = ymax - (y_center_km * 1.e3)
    dists = np.sqrt((horiz_prof[:, x_col] - x_center) ** 2 + (horiz_prof[:, y_col] - y_center_model) ** 2)
    i_near = np.argmin(dists)

    vx_near = horiz_prof[i_near, vx_col]
    vy_near = horiz_prof[i_near, vy_col]
    dip_rad = np.deg2rad(dip_deg)
    vs_near = vx_near * np.cos(dip_rad) - vy_near * np.sin(dip_rad)

    return vs_near, vx_near, vy_near


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
saved_stresses_name = ''.join(['text_files/TESTC/',model_name,'.z',str(analysis_depth/1.e3),'.shear-dz',str(analysis_depth_dz/1.e3),'.ds',str(ds/1.e3),'.prof-dz',str(dz/1.e3),'km.txt'])
os.makedirs(os.path.dirname(saved_stresses_name), exist_ok=True)
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

# create low res grid for plotting
xmin_plot=0; ymin_plot=500.e3; grid_res=2.5e3
X_low, Y_low = create_grid(xmin_plot,xmax,ymin_plot,ymax,grid_res)

# higher res grid for detailed calculations
xmin_plot2=2500.e3; xmax_plot2=3750.e3
ymin_plot2=ymax-600.e3; grid_res2=0.25e3
X_low2, Y_low2 = create_grid(xmin_plot2,xmax_plot2,ymin_plot2,ymax,grid_res2)

first_time=8
saved_stresses = np.zeros(((max_time-first_time),30)) 
                                                     
ind = 0 

for time in range(first_time,max_time,1):

    # plot name
    csv_filename=''.join([csvs_loc,model_name,'/full.',str(time),'.csv'])
    plotname=''.join([plot_loc,'/',str(time),'.png'])
    curvature_plotname=''.join([plot_loc,'/curvature.',str(time),'.png'])

    print("-------")
    print("t = %.0f" % (time))
    model_data  = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

    # print "interpolating model outputs to regular grid..."
    visc   = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,visc_col], (X_low, Y_low), method='linear')
    llith  = griddata((model_data[:,x_col], model_data[:,y_col]), model_data[:,c_llith_col], (X_low2, Y_low2), method='linear')

    # get lithosphere contour, and then trim it to get slab mid-plane
    llith_cont = plt.contour(X_low2/1.e3, (ymax-Y_low2)/1.e3, llith, levels=[0.5])
    print(llith_cont)
    llith_points_tmp = llith_cont.allsegs[0][0]
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

    # near surface prof for plate velocities
    analysis_depth_surf      = 30e3 
    surf_prof                = extract_horiz_prof(model_data,analysis_depth_surf,ymax,0.5e3,x_col,y_col)
    print("got surf prof")
    vc, vsp, vop, trench_loc = get_platevels_from_horiz_prof(surf_prof,x_col,y_col,c_crust_col,vx_col,ymax)
    print("vc = %.1f cm/yr; vsp = %.1f cm/yr; vop = %.1f cm/yr; trench at %.0f km" % (vc,vsp,vop,trench_loc/1e3)) 

    # get dip of mid-plane (and smooth)
    dips_unsmoothed = get_dip_slab_midplane(llith_points) # degrees
    dips = savgol_filter(dips_unsmoothed[:,0],601,3)
    K, dK, K_unsmoothed, dK_unsmoothed = get_curvature_slab_midplane(llith_points,dips) # K [rads/m], dK [rads/m^2]

    # get dips at the relevant depths
    dip_midmant,       xloc_dip,       sloc_dip 		= get_dip_at_certain_depth(dips,llith_points,analysis_depth)
    dip_midmant_shall, xloc_dip_shall, sloc_dip_shall 	= get_dip_at_certain_depth(dips,llith_points,analysis_depth_shall)
    dip_midmant_deep,  xloc_dip_deep,  sloc_dip_deep  	= get_dip_at_certain_depth(dips,llith_points,analysis_depth_deep)

    # get curvature and dcurvature at relevant depth
    K_midmant, dK_midmant, xloc_K, sloc_K                           = get_curvature_at_certain_depth(K,dK,llith_points,analysis_depth)
    K_midmant_shall, dK_midmant_shall, xloc_K_shall, sloc_K_shall   = get_curvature_at_certain_depth(K,dK,llith_points,analysis_depth_shall)
    K_midmant_deep,  dK_midmant_deep,  xloc_K_deep,  sloc_K_deep    = get_curvature_at_certain_depth(K,dK,llith_points,analysis_depth_deep)

    # slab normal thicknesses at the relevant depths
    slabnorm_thick       = (x_right-x_left)             * np.sin(np.deg2rad(dip_midmant))           # m
    slabnorm_thick_ll    = (x_center-x_left)             * np.sin(np.deg2rad(dip_midmant))      
    slabnorm_thick_ul    = (x_right-x_center)             * np.sin(np.deg2rad(dip_midmant))       
    # shallow
    slabnorm_thick_shall    = (x_right_shall-x_left_shall) * np.sin(np.deg2rad(dip_midmant_shall))  # m
    slabnorm_thick_shall_ll = (x_center_shall-x_left_shall) * np.sin(np.deg2rad(dip_midmant_shall)) 
    slabnorm_thick_shall_ul = (x_right_shall-x_center_shall) * np.sin(np.deg2rad(dip_midmant_shall)) 
    # deep
    slabnorm_thick_deep     = (x_right_deep -x_left_deep)  * np.sin(np.deg2rad(dip_midmant_deep))   # m
    slabnorm_thick_deep_ll  = (x_center_deep -x_left_deep)  * np.sin(np.deg2rad(dip_midmant_deep)) 
    slabnorm_thick_deep_ul  = (x_right_deep -x_center_deep)  * np.sin(np.deg2rad(dip_midmant_deep)) 

    # compute theoretical DP
    DP_anal = drho * 9.81 * slabnorm_thick * np.cos(np.deg2rad(dip_midmant)) # Pa

    # pull out model DP (and deviatoric normal stress) across slab
    Pleft, Pright, Sxxleft, Sxxright, Sxyleft, Sxyright, Pleft_x, Pright_x, Pleft_y, Pright_y, slab_x_left, slab_y_left, slab_x_right, slab_y_right = \
        get_nearslab_stresses(y_center,slabnorm_thick_ul,slabnorm_thick_ll,dip_midmant,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz) # Pa

    Pref = get_farfieldP(y_center,ymax,dz,model_data,P_col,x_col,y_col)
    print("P ref = %.1f" % (Pref/1.e6))

    DP_mod = Pleft - Pright # Pa
    Pleft = Pleft - Pref
    Pright = Pright - Pref

    # DP_missing = DP_anal - DP_mod
    print("DPs: analytical = %.1f, model %.1f" % (DP_anal/1.e6,DP_mod/1.e6))
     
    # other depths for shear stress
    Pleft_shall, Pright_shall, Sxxleft_shall, Sxxright_shall, Sxyleft_shall, Sxyright_shall, Pleft_x_shall, Pright_x_shall, Pleft_y_shall, Pright_y_shall, slab_x_left_shall, slab_y_left_shall, slab_x_right_shall, slab_y_right_shall = \
        get_nearslab_stresses(y_center_shall,slabnorm_thick_ul,slabnorm_thick_ll,dip_midmant_shall,model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz) # Pa
    Pleft_deep, Pright_deep, Sxxleft_deep, Sxxright_deep, Sxyleft_deep, Sxyright_deep, Pleft_x_deep, Pright_x_deep, Pleft_y_deep, Pright_y_deep, slab_x_left_deep, slab_y_left_deep, slab_x_right_deep, slab_y_right_deep = \
        get_nearslab_stresses(y_center_deep, slabnorm_thick_ul,slabnorm_thick_ll, dip_midmant_deep, model_data,ds,ymax,x_col,y_col,P_col,c_llith_col,c_ulith_col,sxx_col,sxy_col,dz) # Pa 
    DP_mod_shall = Pleft_shall - Pright_shall;  DP_mod_deep  = Pleft_deep  - Pright_deep

    print("got near slab stresses...")

    # get slab stresses along slab-perp profiles
    sxx_slab, sxy_slab, P_slab, profile = \
        get_stress_profile(model_data,slab_x_left,slab_y_left,slab_x_right,slab_y_right,ymax,x_col,y_col,sxx_col,sxy_col,P_col)
    sxx_slab_shall, sxy_slab_shall, P_slab_shall, profile_shall = \
        get_stress_profile(model_data,slab_x_left_shall,slab_y_left_shall,slab_x_right_shall,slab_y_right_shall,ymax,x_col,y_col,sxx_col,sxy_col,P_col)
    sxx_slab_deep, sxy_slab_deep, P_slab_deep, profile_deep = \
        get_stress_profile(model_data,slab_x_left_deep,slab_y_left_deep,slab_x_right_deep,slab_y_right_deep,ymax,x_col,y_col,sxx_col,sxy_col,P_col)

    # rotate to get shear stresses in slab
    slab_shear_stress       = convert_to_slabnorm_shearstress(sxx_slab,      sxy_slab,      dip_midmant)            # sxx = "shear stress:0" and syy = "shear stress:4" (+VE = compression).  sxy = "shear stress:1" 
    slab_shear_stress_shall = convert_to_slabnorm_shearstress(sxx_slab_shall,sxy_slab_shall,dip_midmant_shall)
    slab_shear_stress_deep  = convert_to_slabnorm_shearstress(sxx_slab_deep, sxy_slab_deep, dip_midmant_deep)

    # rotate to get slab-normal normal stress in slab
    slab_norm_devstress,       slab_norm_fullstress         = convert_to_slabnorm_normstress(sxx_slab, sxy_slab, P_slab,  dip_midmant)
    slab_norm_devstressB,      slab_norm_fullstressB        = convert_to_slabnorm_normstress_TEST(sxx_slab, sxy_slab, P_slab,  dip_midmant)
   

    # get deviatoric normal stress above and below slab (to potentially add to DP)
    Snorm_left  = Sxyleft*np.sin(np.deg2rad(2.*dip_midmant))  + Sxxleft*((np.sin(np.deg2rad(dip_midmant)))**2 - (np.cos(np.deg2rad(dip_midmant)))**2)
    Snorm_right = Sxyright*np.sin(np.deg2rad(2.*dip_midmant)) + Sxxright*((np.sin(np.deg2rad(dip_midmant)))**2 - (np.cos(np.deg2rad(dip_midmant)))**2)
    Snorm_contrib = Snorm_left-Snorm_right # positive = upwards push (helping to support buoyancy) 
    print("Dev stresses at the slab surfaces: %.3f MPa upwards push" % (Snorm_contrib/1.e6))

    # cut out the boundaries of the slabs (which have stress discontinuities)
    profcut_min = 10
    profcut_max = 10
    inds                            = np.where((profile[:,2]>profcut_min) & (profile[:,2]<profile[len(profile)-1,2]-profcut_max))[0]
    slab_shear_stress_cut           = slab_shear_stress[inds]
    slab_norm_devstress_cut         = slab_norm_devstress[inds]
    slab_norm_devstress_cutB        = slab_norm_devstressB[inds]
    profile_cut                     = profile[inds,:]
    inds_shall                      = np.where((profile_shall[:,2]>profcut_min) & (profile_shall[:,2]<profile_shall[len(profile_shall)-1,2]-profcut_max))[0]
    slab_shear_stress_shall_cut     = slab_shear_stress_shall[inds_shall]
    profile_shall_cut           = profile_shall[inds_shall,:]
    inds_deep                   = np.where((profile_deep[:,2]>profcut_min) & (profile_deep[:,2]<profile_deep[len(profile_deep)-1,2]-profcut_max))[0]
    slab_shear_stress_deep_cut  = slab_shear_stress_deep[inds_deep]
    profile_deep_cut            = profile_deep[inds_deep,:]

    # integrate stresses to get slab-perp shear force
    slab_force          = trapz(slab_shear_stress_cut[:,0],         profile_cut[:,2]*1.e3)
    slab_force_shall    = trapz(slab_shear_stress_shall_cut[:,0],   profile_shall_cut[:,2]*1.e3)
    slab_force_deep     = trapz(slab_shear_stress_deep_cut[:,0],    profile_deep_cut[:,2]*1.e3)

    # get ds 
    ds_shall = np.sqrt((x_center - x_center_shall)**2 + (y_center*1e3 - y_center_shall*1e3)**2)
    ds_deep  = np.sqrt((x_center - x_center_deep)**2  + (y_center*1e3 - y_center_deep*1e3)**2)

    # get d(shear stress)/ds 
    slab_stress_term = (slab_force_deep - slab_force_shall)/(ds_shall + ds_deep)
    slab_stress_term_b = (slab_force_deep - slab_force)/(ds_deep)
    slab_stress_term_c = (slab_force - slab_force_shall)/(ds_shall)

    print("slab shear stress term = %.2f MPa (others: %.2f, %.2f)" % (slab_stress_term/1.e6,slab_stress_term_b/1.e6,slab_stress_term_c/1.e6))

    # velocity- and curvature-length scales for dQ/ds scaling:
    # dQ/ds ~ eta H K v_s (2/L_v + 1/L_K), with
    # L_v = |v_s / (dv_s/ds)| and L_K = |K / (dK/ds)|
    vs_mid, vx_mid, vy_mid = get_vs_near_slab_center(midmant_prof, x_center, y_center, ymax, x_col, y_col, vx_col, vy_col, dip_midmant)
    vs_shall, vx_shall, vy_shall = get_vs_near_slab_center(midmant_prof_shall, x_center_shall, y_center_shall, ymax, x_col, y_col, vx_col, vy_col, dip_midmant_shall)
    vs_deep, vx_deep, vy_deep = get_vs_near_slab_center(midmant_prof_deep, x_center_deep, y_center_deep, ymax, x_col, y_col, vx_col, vy_col, dip_midmant_deep)

    dvs_ds = (vs_deep - vs_shall) / (ds_shall + ds_deep)

    eps = 1.e-30
    if np.abs(dvs_ds) > eps and np.abs(vs_mid) > eps:
        Lv = np.abs(vs_mid / dvs_ds)
    else:
        Lv = np.nan

    if np.abs(dK_midmant) > eps and np.abs(K_midmant) > eps:
        Lk = np.abs(K_midmant / dK_midmant)
    else:
        Lk = np.nan

    slab_visc_mid, slab_visc_x_km, slab_visc_y_km = get_slabvisc_from_horiz_prof(midmant_prof, x_col, y_col, visc_col, x_center, y_center, ymax)

    if np.isfinite(Lv) and np.isfinite(Lk) and Lv > 0 and Lk > 0:
        L_total_inv = (2.0 / Lv) + (1.0 / Lk)
        dQds_scaling_splitL = slab_visc_mid * slabnorm_thick * (K_midmant * vs_mid) * L_total_inv
    else:
        L_total_inv = np.nan
        dQds_scaling_splitL = np.nan

    print("L_v = %.1f km, L_K = %.1f km, scaled dQ/ds = %.2f MPa" % (Lv/1.e3, Lk/1.e3, dQds_scaling_splitL/1.e6))

    # integrate normal stresses through the slab
    slabnorm_stress                 = trapz(slab_norm_devstress_cut[:,0],       profile_cut[:,2]*1.e3) 
    slabnorm_stress_fullterm        = slabnorm_stress * K_midmant
    slabnorm_stressB                = trapz(slab_norm_devstress_cutB[:,0],       profile_cut[:,2]*1.e3) 
    slabnorm_stress_fullterm_TEST   = slabnorm_stressB * K_midmant   
    print("slab norm stress = %.2f MPa m, slab norm stress term = %.2f MPa" % (slabnorm_stress/1.e6,slabnorm_stress_fullterm/1.e6))

    # # calculate curvature term for dp
    # curvature_term = 0.5*(slabnorm_thick*K_midmant) * (Pleft+Pright)

    saved_stresses[ind,:] = time, DP_mod_shall, DP_mod_deep, DP_mod, DP_anal, dip_midmant, slab_stress_term, slab_stress_term_b,  \
                            slab_stress_term_c,  slabnorm_thick, Snorm_contrib, K_midmant, dK_midmant, K_midmant_shall, K_midmant_deep, \
                            Pleft, Pright, slabnorm_stress_fullterm, slabnorm_stress_fullterm_TEST, vc, vsp, \
                            vs_mid, dvs_ds, Lv, Lk, dQds_scaling_splitL, slab_visc_mid, slab_visc_x_km, slab_visc_y_km, L_total_inv
    ind = ind + 1

    ###################### plotting - 1 #########################
    fig=plt.figure()
    gs=GridSpec(6,1) 

    # plot viscosity field
    ax1=fig.add_subplot(gs[0,0])
    visc_plot = ax1.contourf(X_low/1.e3, (ymax-Y_low)/1.e3, np.log10(visc), cmap=cm.get_cmap('plasma_r'),levels=np.linspace(19,24,501))
    ax1.set_ylim([(ymax-ymin_plot)/1.e3,0])   
    ax1.tick_params(direction='out',length=2, labelsize=6)
    ax1.annotate(''.join(['DP: mod = ',str("%.1f" % (DP_mod/1.e6)),', anal = ',str("%.1f" % (DP_anal/1.e6)),' MPa, vc = ',str("%.1f" % (vc)),' cm/yr']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=6.5,color='k')        
    ax1.plot(llith_points[:,0], llith_points[:,1], linewidth=0.5, color='black',zorder=6)
    ax1.plot(profile_shall[:,0],profile_shall[:,1], linewidth=0.5, color='green',zorder=6)
    ax1.plot(profile_deep[:,0], profile_deep[:,1],  linewidth=0.5, color='blue',zorder=6)
    ax1.plot(profile[:,0],      profile[:,1],       linewidth=0.5, color='red',zorder=6)
    ax1.scatter(Pleft_x/1.e3,  Pleft_y/1.e3,  s=0.25,color='black',zorder=3)
    ax1.scatter(Pright_x/1.e3, Pright_y/1.e3, s=0.25,color='black',zorder=3)
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

    # plot in-slab shear stress
    ax5=fig.add_subplot(gs[4,0])
    ax5.set_ylabel(r"$\tau_{s,n}$ [MPa]",size=6.5)
    ax5.set_xlim(-10,100); 
    ax5.tick_params(axis='x', labelsize=6)
    ax5.tick_params(axis='y', labelsize=6)
    ax5.plot(profile[:,2],      slab_shear_stress/1.e6,      color='red', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
    ax5.plot(profile_shall[:,2],slab_shear_stress_shall/1.e6, color='green', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
    ax5.plot(profile_deep[:,2], slab_shear_stress_deep/1.e6,  color='blue', linestyle='-', linewidth=1,zorder=2,alpha=0.4)  
    ax5.plot(profile_cut[:,2],      slab_shear_stress_cut/1.e6,       color='red', linestyle='-', linewidth=1.5,zorder=3)  
    ax5.plot(profile_shall_cut[:,2],slab_shear_stress_shall_cut/1.e6, color='green', linestyle='-', linewidth=1.5,zorder=3)  
    ax5.plot(profile_deep_cut[:,2], slab_shear_stress_deep_cut/1.e6,  color='blue', linestyle='-', linewidth=1.5,zorder=3) 
    misfit = DP_mod-DP_anal-slab_stress_term-slabnorm_stress_fullterm
    ax5.annotate(''.join(['slab shear stress = ',str("%.1f" % (slab_stress_term/1.e6)),', misfit = ',str("%.1f" % (misfit/1.e6)),' MPa']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=5,color='k')        
    ax5.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)

    # plot in-slab norm stress
    ax6=fig.add_subplot(gs[5,0])
    ax6.set_ylabel(r"$\sigma_{s,s}$ [MPa]",size=6.5)
    ax6.set_xlim(-10,100); 
    ax6.tick_params(axis='x', labelsize=6)
    ax6.tick_params(axis='y', labelsize=6)
    ax6.plot(profile[:,2],slab_norm_devstress/1.e6, color='gray', linestyle='-', linewidth=1.5,zorder=3)  
    ax6.plot(profile_cut[:,2],slab_norm_devstress_cut/1.e6,       color='black', linestyle='-', linewidth=1.5,zorder=3)  
    misfit = DP_mod-DP_anal-slab_stress_term-slabnorm_stress_fullterm
    ax6.axhline(y=0, color='gray',linestyle='--',linewidth=1, zorder=1)
    ax6.annotate(''.join(['slab norm stress term = ',str("%.1f" % (slabnorm_stress_fullterm/1.e6)), 'MPa, norm stress = ',str("%.1f" % (slabnorm_stress/1.e6)),'MPa m, K = ',str("%.5f" % (K_midmant*1000)),' 1/km)']), xy=(0.01,0.12), xycoords='axes fraction',verticalalignment='center',horizontalalignment='left',fontsize=5,color='k')        

    # print("saving fields figure to %s..." % plotname)
    plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
    plt.clf()


np.savetxt(saved_stresses_name, saved_stresses)
