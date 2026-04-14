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
from functions import make_strictly_ascending
from functions import load_data_file
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



plotname=''.join(['plots/maps_dip-age.png'])
plotname_pdf=''.join(['plots/maps_dip-age.pdf'])

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
vminDip, vmaxDip = 25, 90
norm1 = mcolors.Normalize(vmin=vminDip, vmax=vmaxDip)
vminAge, vmaxAge = 0, 200
norm2 = mcolors.Normalize(vmin=vminAge, vmax=vmaxAge)
DP_array = []
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

    if not np.isnan(dip_shall) and age < 250 and not np.isnan(vc) and not np.isnan(K):

        x1, y1 = m1(lon_center, lat_center)
        ck = ax1.scatter(x1, y1, s=23, c=dip_shall, cmap='BrBG',      norm=norm1, edgecolors='black', linewidths=0.4, zorder=10)
        cv = ax2.scatter(x1, y1, s=23, c=age,       cmap='inferno_r', norm=norm2, edgecolors='black', linewidths=0.4, zorder=10)
                
        # store DP and stress scaling
        DP_array.append((0, np.abs(K)))

DP_array = np.array(DP_array)

ck_bar = plt.colorbar(ck, ax=ax1, label=r'$\theta$  [$\degree$]',  shrink=0.5, pad=0.05)
ck_bar.set_ticks([30, 40, 50, 60, 70, 80, 90])
cv_bar = plt.colorbar(cv, ax=ax2, label=r'age  [Ma]',  shrink=0.5, pad=0.05)

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


plt.tight_layout()
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=700)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()
