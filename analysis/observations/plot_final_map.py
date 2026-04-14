#!/usr/bin/env python
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as mcolors
import numpy as np
from netCDF4 import Dataset
import sys
import glob
import math
from scipy.interpolate import RegularGridInterpolator
from functions import read_pb2002_boundaries, plot_ocean_age
from functions import haversine, calculate_bearing, destination_point
from functions import make_strictly_ascending, compute_DP_hs, compute_DP_pl, compute_H_eff
from functions import load_data_file, stats_data_file, stats_DP
import matplotlib.font_manager as fm
font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
myriad_pro = fm.FontProperties(fname=font_path)


mpl.rcParams['font.family'] = 'Myriad Pro'  # Now it should work if properly installed!
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.labelsize'] = 10
mpl.rcParams['axes.labelpad'] = 1.25
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.pad'] = 2
mpl.rcParams['ytick.major.pad'] = 2
mpl.rcParams['xtick.major.size'] = 2.5
mpl.rcParams['ytick.major.size'] = 2.5
mpl.rcParams['xtick.minor.size'] = 1.25
mpl.rcParams['ytick.minor.size'] = 1.25
scaling_thresh = 5  # MPa, threshold for scaling parameter


# constants
slab_visc = float(sys.argv[1])    # Pa s [ref=4e22]
alpha = float(sys.argv[2])        # 1/K  [ref=3.28e-5]
Tm = float(sys.argv[3])           # degC [ref=1333]
diffusivity = float(sys.argv[4])  # m^2/s [ref=8.044e-7]
plate_thick = float(sys.argv[5])  # m [ref=88e3]
crust_thick = float(7e3)          # m [ref=7e3]
cooling_model = 'plate-cooling'

# conversion factors
cmyr_to_ms = 1e-2 / 3.154e7
Ma_to_s = 1e6 * 3.154e7


plotname=''.join(['plots/maps.slab',str(slab_visc),'.alpha',str(alpha),'.T',str(Tm),'.k',str(diffusivity),'.platethick',str(plate_thick),'.',cooling_model,'.png'])
plotname_pdf=''.join(['plots/maps.slab',str(slab_visc),'.alpha',str(alpha),'.T',str(Tm),'.k',str(diffusivity),'.platethick',str(plate_thick),'.',cooling_model,'.pdf'])

# ------------------------------------------------
# --- Set up the map ---
# ------------------------------------------------
fig = plt.figure(figsize=(8,6))
G = gridspec.GridSpec(2,3)

# First subplot
ax1 = plt.subplot(G[0, 0:2])
m1 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax1)
m1.drawcoastlines(linewidth=0)
m1.drawmapboundary(fill_color='white')
m1.fillcontinents(color='silver', lake_color='silver')

# Second subplot
ax2 = plt.subplot(G[1, 0:2])
m2 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax2)
m2.drawcoastlines(linewidth=0)
m2.drawmapboundary(fill_color='white')
m2.fillcontinents(color='silver', lake_color='silver')

#------------------------------------------------
# ------ plot slab contours ---------------------
#------------------------------------------------
slab_files = glob.glob("data/Slab2/*_dep_*grd")
for grd_file in slab_files:
    nc = Dataset(grd_file, 'r')
    lon_arr = nc.variables['x'][:]  # Longitudespython3 plot_final_map.py 4e22 2e-5 1300 1e-6

    lat_arr = nc.variables['y'][:]  # Latitudes
    slab_data = nc.variables['z'][:]  # Slab depth data

    # Make the arrays strictly ascending (in case of duplicate or constant values)
    lat_arr = make_strictly_ascending(lat_arr)
    lon_arr = make_strictly_ascending(lon_arr)

    # close the file
    nc.close()

    interp = RegularGridInterpolator(
        (lat_arr, lon_arr), slab_data, bounds_error=False, fill_value=np.nan)
    
    # Plot the slab data
    lon2d, lat2d = np.meshgrid(lon_arr, lat_arr)
    x2d, y2d = m1(lon2d, lat2d)
    levels = np.arange(-600, 1, 100) 
    cs = ax1.contour(x2d, y2d, slab_data, levels=levels, colors='red', linestyles='solid', linewidths=0.5,zorder=10)

    x2d2, y2d2 = m2(lon2d, lat2d)
    ax2.contour(x2d2, y2d2, slab_data, levels=levels, colors='red', linestyles='solid', linewidths=0.5,zorder=10)

# ------------------------------------------------
# -------- Read in the segment data --------------
# ------------------------------------------------
# Load data from file
data = np.genfromtxt("data/segment_data.txt", delimiter=',', dtype=str)
data[data == 'None'] = np.nan #  lon_center, lat_center, shallow_K, deep_K, chosen_vc, age, lon_for_age, lat_out, dip, dip_shall
segment_data = np.array(data, dtype=float)

# ------------------------------------------------
# ---- compute and plot DP and K vc eta ----------
# ------------------------------------------------
DPmin, DPmax = 10, 60
vminK, vmaxK = 0, 15
norm1 = mcolors.Normalize(vmin=DPmin, vmax=DPmax)
norm2 = mcolors.Normalize(vmin=vminK, vmax=vmaxK)
DP_array = []
min_dp=1e9
max_dp=0
for i in range(len(segment_data)):

    # load everything in
    lon_center, lat_center = segment_data[i,0], segment_data[i,1]
    age = segment_data[i,5]               # Ma
    dip_deep  = segment_data[i, 8]        # degrees
    dip_shall = segment_data[i, 9]        # degrees
    shallow_K = segment_data[i, 2]        # 1/km
    deep_K = segment_data[i, 3]           # 1/km
    vc = segment_data[i,4]                # cm/yr
    if np.isnan(shallow_K) and np.isnan(deep_K):
        K = np.nan
    elif np.isnan(deep_K):
        K = shallow_K/1000.   # now 1/m
    else:
        K = deep_K/1000.     

    x1, y1 = m1(lon_center, lat_center)
    x2, y2 = m2(lon_center, lat_center)

    if not np.isnan(dip_shall) and age < 250 and not np.isnan(vc) and not np.isnan(K):

        # compute/plot ηHKvc/L_eff -------------------
        H_eff = compute_H_eff(age, Tm, k=diffusivity, plate_thick=plate_thick)  # m
        stress_scaling = K * (vc * cmyr_to_ms) * slab_visc * H_eff / 1.624e6 * 1e-6  # MPa
        x2, y2 = m2(lon_center, lat_center)
        if np.abs(stress_scaling) > scaling_thresh:
            edgecolor = 'gray'
            edgethick = 0.3
            zord = 10
        else:
            edgecolor = 'black'
            edgethick = 0.35
            zord = 11
        ck = ax2.scatter(x2, y2, s=23, c=np.abs(stress_scaling), cmap='BrBG', norm=norm2, edgecolors=edgecolor, linewidths=edgethick, zorder=zord)

        # compute/plot DP ------------------------
        if not np.isnan(dip_deep):
            DP = compute_DP_pl(age, dip_deep, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
            print(lon_center,lat_center,age,dip_deep,DP)
        else:
            DP = compute_DP_pl(age, dip_shall, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
                
        if DP < min_dp:
            min_dp=DP
        if DP > max_dp:
            max_dp=DP
        cp = ax1.scatter(x1, y1, s=23, c=DP, cmap='plasma_r', norm=norm1, edgecolors=edgecolor, linewidths=edgethick, zorder=zord)

        # store DP and stress scaling
        DP_array.append((DP, np.abs(stress_scaling)))

print(min_dp,max_dp)

DP_array = np.array(DP_array)

cp_bar = plt.colorbar(cp, ax=ax1, label=r'$\Delta P$   [MPa]', extend='max', shrink=0.5, pad=0.05)
cp_bar.set_ticks([10, 20, 30, 40, 50, 60])
ck_bar = plt.colorbar(ck, ax=ax2, label=r'($\eta H K V_{C}$)/$L_\mathrm{eff}$  [MPa]', extend='max', shrink=0.5, pad=0.05)

# ------------------------------------------------
# --- Plot plate boundaries ---
# ------------------------------------------------
boundaries = read_pb2002_boundaries("data/PB2002/PB2002_boundaries.dig.txt")
for seg in boundaries:
    if len(seg['coords']) < 2:
        continue  # need at least two points to plot a line
    lons, lats = zip(*seg['coords'])
    x, y = m1(lons, lats)
    if seg['symbol'] in ["/", "\\"]:
        color = 'magenta'
        lw = 1.25
    else:
        color = 'black'
        lw = 0.5
    # Handle boundaries that cross the prime meridian.
    if min(lons) < 0 and max(lons) > 0:
        lons1, lats1, lons2, lats2 = [], [], [], []
        for i in range(len(lons)):
            if lons[i] > 0:
                lons1.append(lons[i])
                lats1.append(lats[i])
            else:
                lons2.append(lons[i])
                lats2.append(lats[i])
        m1.plot(*m1(lons1, lats1), color=color, linewidth=lw)
        m1.plot(*m1(lons2, lats2), color=color, linewidth=lw)
        m2.plot(*m2(lons1, lats1), color=color, linewidth=lw)
        m2.plot(*m2(lons2, lats2), color=color, linewidth=lw)
    else:
        m1.plot(x, y, color=color, linewidth=lw)
        m2.plot(x, y, color=color, linewidth=lw)

# ------------------------------------------------
# ----------- Plot DP histogram ------------------
# ------------------------------------------------

# Define the independent variables:
T_vals = np.arange(1100, 1450+10, 10)
plate_thick_vals = np.arange(75e3, 135e3+2e3, 2e3)

# Initialize a grid to store the DPmean_thresh values.
# The grid shape will be (# of plate_thick_vals, # of T_vals)
DPmean_thresh_grid = np.zeros((len(plate_thick_vals), len(T_vals)))

# Loop over all combinations of T and plate thickness and compute DPmean_thresh.
for i, plate_thick_val in enumerate(plate_thick_vals):
    for j, T_val in enumerate(T_vals):
        data = load_data_file(slab_visc, alpha, float(T_val), diffusivity, float(plate_thick_val), crust_thick, cooling_model)
        DPmean_thresh, DPmean_tot = stats_DP(data,scaling_thresh)  # threshold of 10 MPa
        DPmean_thresh_grid[i, j] = DPmean_thresh

        
# get reference mean
#data_ref = load_data_file(slab_visc, alpha, Tm,  diffusivity, plate_thick, crust_thick, cooling_model)
#DPmean_thresh_ref, DPmean_tot_ref = stats_DP(data_ref,scaling_thresh)  # threshold of 10 MPa
#print(DPmean_thresh_ref, DPmean_tot_ref)
# --

plate_thick_vals_km = plate_thick_vals / 1e3

# Create the subplot (ax3) for the colored grid.
ax3 = plt.subplot(G[0, 2])

contour_levels1 = np.arange(20, 40+2, 2)
contours1 = ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid, levels=contour_levels1, colors='black', linewidths=1, linestyles='solid')
ax3.clabel(contours1, inline=True, fontsize=10)

# plot paramters of the map
ax3.plot(Tm, plate_thick/1e3, 'k*', markersize=9)

# Set axis labels and ticks.
ax3.set_xlabel('Basal temperature  [°C]')
ax3.set_ylabel('Plate thickness  [km]')
ax3.set_xticks(np.arange(1100, 1400+100, 100))
ax3.set_yticks(np.arange(75,135,10))
ax3.tick_params(axis='both', labelsize=10)
ax3.labelsize = 12


# ------------------------------------------------
# -------- Plot scaling histogram ----------------
# ------------------------------------------------

# load in scaling & DP values for other slab viscs
array1  = load_data_file(4e20,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array2  = load_data_file(8e20,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array3  = load_data_file(1e21,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array4  = load_data_file(4e21,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array5  = load_data_file(8e21,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array6  = load_data_file(1e22,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array7  = load_data_file(4e22,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array7b  = load_data_file(6e22,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array8  = load_data_file(8e22,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array9  = load_data_file(1e23,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array10 = load_data_file(4e23,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array11 = load_data_file(8e23,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array12 = load_data_file(1e24,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)
array13 = load_data_file(4e24,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model)

mean1, std1, perc1 = stats_data_file(array1,scaling_thresh)
mean2, std2, perc2 = stats_data_file(array2,scaling_thresh)
mean3, std3, perc3 = stats_data_file(array3,scaling_thresh)
mean4, std4, perc4 = stats_data_file(array4,scaling_thresh)
mean5, std5, perc5 = stats_data_file(array5,scaling_thresh)
mean6, std6, perc6 = stats_data_file(array6,scaling_thresh)
mean7, std7, perc7 = stats_data_file(array7,scaling_thresh)
mean7b, std7b, perc7b = stats_data_file(array7b,scaling_thresh)
mean8, std8, perc8 = stats_data_file(array8,scaling_thresh)
mean9, std9, perc9 = stats_data_file(array9,scaling_thresh)
mean10, std10, perc10 = stats_data_file(array10,scaling_thresh)
mean11, std11, perc11 = stats_data_file(array11,scaling_thresh)
mean12, std12, perc12 = stats_data_file(array12,scaling_thresh)
mean13, std13, perc13 = stats_data_file(array13,scaling_thresh)

ax4 = plt.subplot(G[1, 2])

color='white'
edgescolor='black'
size=70
ax4.scatter(np.log10(4e20), perc1,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(8e20), perc2,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(1e21), perc3,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(4e21), perc4,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(8e21), perc5,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(1e22), perc6,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(4e22), perc7,  s=size, c='black', edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(6e22), perc7b,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(8e22), perc8,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(1e23), perc9,  s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(4e23), perc10, s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(8e23), perc11, s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(1e24), perc12, s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)
ax4.scatter(np.log10(4e24), perc13, s=size, c=color, edgecolors=edgescolor, linewidths=0.65, zorder=10)

# ticks on bottom
ax4.set_xticks([20, 21, 22, 23, 24, 25])
ax4.set_xticklabels(['$10^{20}$', '$10^{21}$', '$10^{22}$', '$10^{23}$', '$10^{24}$', '$10^{25}$'])
ax4.xaxis.set_ticks_position('bottom')

# ticks on top
ax_top = ax4.twiny()  # Create a new x-axis on the top
ax_top.set_xlim(ax4.get_xlim())
top_ticks = [np.log10(2e20), np.log10(2e21), np.log10(2e22), np.log10(2e23), np.log10(2e24)]
top_tick_labels = ['1', '10', '100', '1000', '10000']
ax_top.set_xticks(top_ticks)
ax_top.set_xticklabels(top_tick_labels)
side_ticks = [0, 20, 40, 60, 80, 100]
side_tick_labels = ['0', '20 %', '40 %', '60 %', '80 %', '100 %']
ax_top.set_yticks(side_ticks)
ax_top.set_yticklabels(side_tick_labels)

ax4.set_ylabel(r'($\eta H K V_{C}$)/$L_\mathrm{eff}$  <  5 MPa')
ax4.set_xlabel(r'$\eta$  [Pa s]' )
ax_top.set_xlabel(r'$\eta$ / $\eta_{mantle}$')
ax4.set_xlim(20, 25)
ax4.set_ylim(-4, 104)

ax4.tick_params(axis='both', labelsize=10)
ax_top.tick_params(axis='both', labelsize=10)
ax4.labelsize = 12

plt.tight_layout()
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=700)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()
