#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as mcolors
import numpy as np
from netCDF4 import Dataset
import glob
import math
from scipy.interpolate import RegularGridInterpolator
from functions import read_pb2002_boundaries, plot_ocean_age
from functions import haversine, calculate_bearing, destination_point
from functions import make_strictly_ascending

# ------------------------------------------------
# --- Set up the map ---
# ------------------------------------------------
fig = plt.figure(figsize=(12,6))

# First subplot
ax1 = fig.add_subplot(111)
m1 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax1)
m1.drawcoastlines(linewidth=0.5)
m1.drawmapboundary(fill_color='lightblue')
m1.fillcontinents(color='white', lake_color='lightblue')


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
        lw = 2
    else:
        color = 'black'
        lw = 1
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
    else:
        m1.plot(x, y, color=color, linewidth=lw)

# ---------------------------------------------------------
# --- Pre-load slab grid data from Slab2 files and plot ---
# ---------------------------------------------------------
# Build a list of slab grid dictionaries for later lookup.

slab_files = glob.glob("data/Slab2/*_dep_*grd")
slab_grids = []
for grd_file in slab_files:
    nc = Dataset(grd_file, 'r')
    lon_arr = nc.variables['x'][:]  # Longitudes
    lat_arr = nc.variables['y'][:]  # Latitudes
    slab_data = nc.variables['z'][:]  # Slab depth data

    # Make the arrays strictly ascending (in case of duplicate or constant values)
    lat_arr = make_strictly_ascending(lat_arr)
    lon_arr = make_strictly_ascending(lon_arr)

    # close the file
    nc.close()

    interp = RegularGridInterpolator(
        (lat_arr, lon_arr), slab_data, bounds_error=False, fill_value=np.nan)
    
    grid_info = {
        'lon_arr': lon_arr,
        'lat_arr': lat_arr,
        'slab_data': slab_data,
        'interp': interp,
        'min_lon': np.min(lon_arr),
        'max_lon': np.max(lon_arr),
        'min_lat': np.min(lat_arr),
        'max_lat': np.max(lat_arr)
    }

    slab_grids.append(grid_info)

    # Plot the slab data
    # Create 2D grids for lon and lat, then project to map coordinates
    lon2d, lat2d = np.meshgrid(lon_arr, lat_arr)
    x2d, y2d = m1(lon2d, lat2d)
    levels = np.linspace(-600, 0, 13)  # 51 boundaries produce 50 contour intervals
    cf = ax1.contourf(x2d, y2d, slab_data, levels=levels, cmap='viridis', extend='both')

plt.colorbar(cf, label='Slab Depth', extend='both')


# plt.colorbar(cf2, label='dip')

# ------------------------------------------------
# segments and ages
# ------------------------------------------------
segment_length = 300  # km for splitting the boundary
project_distance = 200  # km to project outboard
project_distance2 = 600  # km to project outboard

# --- Load ocean age dataset for interpolation ---
age_nc = Dataset("data/Muller2008/age.3.6.nc", 'r')
lon_age = age_nc.variables['x'][:]
lat_age = age_nc.variables['y'][:]
age_data = age_nc.variables['z'][:]
interp_age = RegularGridInterpolator((lat_age, lon_age), age_data,
                                       bounds_error=False, fill_value=np.nan)
age_nc.close()


# Lists to store projected point map coordinates and their ocean age.
segment_locs = []
segment_dips = []
segment_data = []


for seg in boundaries:
    if seg['symbol'] not in ["/", "\\"]:
        continue  # Only process subduction zones.
    coords = seg['coords']
    if len(coords) < 2:
        continue

    # Compute cumulative distance along the subduction zone.
    cumdist = [0]
    for i in range(1, len(coords)):
        lon1, lat1 = coords[i-1]
        lon2, lat2 = coords[i]
        d = haversine(lon1, lat1, lon2, lat2)
        cumdist.append(cumdist[-1] + d)
    total_length = cumdist[-1]
    num_segments = int(total_length // segment_length)

    for i in range(num_segments):
        # Determine the center of each 200 km segment.
        center_distance = i * segment_length + segment_length / 2.051
        for j in range(len(cumdist) - 1):
            if cumdist[j] <= center_distance <= cumdist[j+1]:
                frac = (center_distance - cumdist[j]) / (cumdist[j+1] - cumdist[j])
                lon_center = coords[j][0] + frac * (coords[j+1][0] - coords[j][0])
                lat_center = coords[j][1] + frac * (coords[j+1][1] - coords[j][1])
                x_center, y_center = m1(lon_center, lat_center)

                # Compute the bearing along the segment.
                bearing = calculate_bearing(coords[j][0], coords[j][1],
                                            coords[j+1][0], coords[j+1][1])
                # Perpendicular bearings.
                if seg['symbol'] == "/":
                    perp = (bearing + 90) % 360
                    slab_perp = (bearing - 90) % 360
                else:
                    perp = (bearing - 90) % 360
                    slab_perp = (bearing + 90) % 360

                # Compute (first) point for extracting ocean age
                lon_out, lat_out = destination_point(lon_center, lat_center, perp, project_distance)
                lon_for_age = lon_out
                if lon_out < lon_age[0]:
                    lon_for_age += 360
                elif lon_out > lon_age[-1]:
                    lon_for_age -= 360    
                age = interp_age((lat_out, lon_for_age))/100 # Convert to Ma

                # try second point for ocean age if needed
                if age > 300:
                    lon_out, lat_out = destination_point(lon_center, lat_center, perp, project_distance2)
                    if lon_out < lon_age[0]:
                        lon_for_age += 360
                    elif lon_out > lon_age[-1]:
                        lon_for_age -= 360    
                    age = interp_age((lat_out, lon_for_age))/100 # Convert to Ma


                # all of the main quantities
                segment_data.append((lon_center, lat_center, 0, 0, 0, age, lon_for_age, lat_out))

# segment_locs = np.array(segment_locs)    
segment_data = np.array(segment_data)

# plot ages and segment centers
norm = mcolors.Normalize(vmin=0, vmax=200)

x_center, y_center = m1(segment_data[:,0], segment_data[:,1])
ax1.plot(x_center, y_center, 'o', color='magenta', markersize=4)
x_out, y_out = m1(segment_data[:,6], segment_data[:,7])
sc = ax1.scatter(x_out, y_out, c=segment_data[:,5], cmap='jet', norm=norm, s=30, edgecolors='k', zorder=5)
cb1 = plt.colorbar(sc, ax=ax1, label='Oceanic Age (Ma)')

# # Plot ocean age points (if desired)
# plot_ocean_age(m1, "data/Muller2008/age.3.6.nc")

plt.tight_layout()
plt.show()
