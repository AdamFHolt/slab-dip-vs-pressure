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


plotname=''.join(['just-maps.slab',str(slab_visc),'.alpha',str(alpha),'.T',str(Tm),'.k',str(diffusivity),'.platethick',str(plate_thick),'.',cooling_model,'.png'])
plotname_pdf=''.join(['just-maps.slab',str(slab_visc),'.alpha',str(alpha),'.T',str(Tm),'.k',str(diffusivity),'.platethick',str(plate_thick),'.',cooling_model,'.pdf'])

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
slab_files = glob.glob("files/*_dep_*grd")
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
data = np.genfromtxt("files/segment_data.txt", delimiter=',', dtype=str)
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
        stress_scaling = K * (vc * cmyr_to_ms) * slab_visc * H_eff / 1.497e6 * 1e-6  # MPa
        x2, y2 = m2(lon_center, lat_center)
        if np.abs(stress_scaling) > scaling_thresh:
            edgecolor = 'gray'
            edgethick = 0
            zord = 10
        else:
            edgecolor = 'black'
            edgethick = 0.35
            zord = 11
        ck = ax2.scatter(x2, y2, s=23, c=np.abs(stress_scaling), cmap='BrBG', norm=norm2, edgecolors=edgecolor, linewidths=edgethick, zorder=zord)

        # compute/plot DP ------------------------
        if not np.isnan(dip_deep):
            DP = compute_DP_pl(age, dip_deep, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
        else:
            DP = compute_DP_pl(age, dip_shall, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
                
        if DP < min_dp:
            min_dp=DP
        if DP > max_dp:
            max_dp=DP
        cp = ax1.scatter(x1, y1, s=23, c=DP, cmap='plasma_r', norm=norm1, edgecolors=edgecolor, linewidths=edgethick, zorder=zord)

        # store DP and stress scaling
        DP_array.append((DP, np.abs(stress_scaling)))


DP_array = np.array(DP_array)

cp_bar = plt.colorbar(cp, ax=ax1, label=r'$\Delta P$   [MPa]', extend='max', shrink=0.5, pad=0.05)
cp_bar.set_ticks([10, 20, 30, 40, 50, 60])
ck_bar = plt.colorbar(ck, ax=ax2, label=r'($\eta H K V_{C}$)/$L_\mathrm{eff}$  [MPa]', extend='max', shrink=0.5, pad=0.05)

# ------------------------------------------------
# --- Plot plate boundaries ---
# ------------------------------------------------
boundaries = read_pb2002_boundaries("files/PB2002_boundaries.dig.txt")
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

plt.tight_layout()
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=700)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()
