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
ax1 = fig.add_subplot(211)
m1 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax1)
m1.drawcoastlines(linewidth=0.5)
m1.drawmapboundary(fill_color='lightblue')
m1.fillcontinents(color='white', lake_color='lightblue')

# Second subplot
ax2 = fig.add_subplot(212)
m2 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax2)
m2.drawcoastlines(linewidth=0.5)
m2.drawmapboundary(fill_color='lightblue')
m2.fillcontinents(color='white', lake_color='lightblue')



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
        m2.plot(*m2(lons1, lats1), color=color, linewidth=lw)
        m2.plot(*m2(lons2, lats2), color=color, linewidth=lw)
    else:
        m1.plot(x, y, color=color, linewidth=lw)
        m2.plot(x, y, color=color, linewidth=lw)

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

# ------------------------------------------------
# now load in dips
# ------------------------------------------------
dip_files = glob.glob("data/Slab2/*_dip_*grd")
dip_grids = []
for grd_file in dip_files:
    nc = Dataset(grd_file, 'r')
    lon_arr = nc.variables['x'][:]  # Longitudes
    lat_arr = nc.variables['y'][:]  # Latitudes
    dip_data = nc.variables['z'][:]  # Dip data (same format as depth grids)

    # Make the arrays strictly ascending
    lat_arr = make_strictly_ascending(lat_arr)
    lon_arr = make_strictly_ascending(lon_arr)
    nc.close()

    interp = RegularGridInterpolator(
        (lat_arr, lon_arr), dip_data, bounds_error=False, fill_value=np.nan)
    
    grid_info = {
        'lon_arr': lon_arr,
        'lat_arr': lat_arr,
        'dip_data': dip_data,
        'interp': interp,
        'min_lon': np.min(lon_arr),
        'max_lon': np.max(lon_arr),
        'min_lat': np.min(lat_arr),
        'max_lat': np.max(lat_arr)
    }
    dip_grids.append(grid_info)

    # Plot the dip data
    # Create 2D grids for lon and lat, then project to map coordinates
    lon2d, lat2d = np.meshgrid(lon_arr, lat_arr)
    x2d, y2d = m2(lon2d, lat2d)
    levels = np.linspace(10, 80)  
    cf2 = ax2.contourf(x2d, y2d, dip_data, levels=levels, cmap='jet')

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

# ------------------------------------------------
# --- Load in convergence rates ------------------ 
# ------------------------------------------------
lallemand_vc = []
with open("data/Lallemand/Lallemand_et_al-2005_G3_dataset.txt", 'r') as f:
    for line in f:
        data_line = line.strip().split()
        lallemand_lon = float(data_line[2]) % 360
        lallemand_vc.append((lallemand_lon,float(data_line[1]),float(data_line[9])))
lallemand_vc = np.array(lallemand_vc)    

# Lists to store projected point map coordinates and their ocean age.
segment_locs = []
segment_dips = []
segment_data = []
target_depth1 = 280
target_depth2 = 320
target_depth3 = 180
target_depth4 = 220

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
        center_distance = i * segment_length + segment_length / 2.0
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

                # Compute point for extracting ocean age
                lon_out, lat_out = destination_point(lon_center, lat_center, perp, project_distance)

                # Interpolate ocean age at both candidate locations.
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


                min_dist1 = 1e9; min_dist2 = 1e9; min_dist3 = 1e9; min_dist4 = 1e9
                slab_depth1 = None; chosen_slab_depth1 = None
                chosen_lon1 = None; chosen_lat1 = None
                slab_depth2 = None; chosen_slab_depth2 = None
                chosen_lon2 = None; chosen_lat2 = None
                slab_depth3 = None; chosen_slab_depth3 = None

                chosen_lon3 = None; chosen_lat3 = None
                slab_depth4 = None; chosen_slab_depth4 = None
                chosen_lon4 = None; chosen_lat4 = None
                dip_val1 = None; dip_val2 = None; 
                dip_val3 = None; dip_val4 = None
                min_dist_threshold = 5 # km

                for dist in np.linspace(0, 1000, 500):

                    lon_search, lat_search = destination_point(lon_center, lat_center, slab_perp, dist)
                    lon_search = lon_search % 360  # Wrap around to 0-360

                    x_search, y_search = m1(lon_search, lat_search)
                    plt.plot(x_search, y_search, 'o', color='red', markersize=2)

                    # Loop over the available slab grids to see if this point falls inside.
                    for grid in slab_grids:
                        # Check if the search point is within the grid's extent.
                        if (grid['min_lon'] <= lon_search <= grid['max_lon'] and
                            grid['min_lat'] <= lat_search <= grid['max_lat']):

                            # Get the slab depth via interpolation.
                            slab_depth  = grid['interp']((lat_search, lon_search))

                            dist_misfit1 = abs((-1.0*target_depth1) - slab_depth)
                            if dist_misfit1 < min_dist1 and dist_misfit1 < min_dist_threshold:
                                min_dist1 = dist_misfit1
                                chosen_slab_depth1 = slab_depth
                                chosen_lon1 = lon_search
                                chosen_lat1 = lat_search
                                chosen_grid1_minlon = np.round(grid['min_lon'])

                            dist_misfit2 = abs((-1.0*target_depth2) - slab_depth)
                            if dist_misfit2 < min_dist2 and dist_misfit2 < min_dist_threshold:
                                min_dist2 = dist_misfit2
                                chosen_slab_depth2 = slab_depth
                                chosen_lon2 = lon_search
                                chosen_lat2 = lat_search
                                chosen_grid2_minlon = np.round(grid['min_lon'])

                            dist_misfit3 = abs((-1.0*target_depth3) - slab_depth)
                            if dist_misfit3 < min_dist3 and dist_misfit3 < min_dist_threshold:
                                min_dist3 = dist_misfit3
                                chosen_slab_depth3 = slab_depth
                                chosen_lon3 = lon_search
                                chosen_lat3 = lat_search
                                chosen_grid3_minlon = np.round(grid['min_lon'])
                            
                            dist_misfit4 = abs((-1.0*target_depth4) - slab_depth)
                            if dist_misfit4 < min_dist4 and dist_misfit4 < min_dist_threshold:
                                min_dist4 = dist_misfit4
                                chosen_slab_depth4 = slab_depth
                                chosen_lon4 = lon_search
                                chosen_lat4 = lat_search
                                chosen_grid4_minlon = np.round(grid['min_lon'])


                # Now, extract dip values at the chosen points from the dip grids.
                dip_val1 = None
                for grid in dip_grids:
                    if (chosen_slab_depth1 != None and grid['min_lon'] <= chosen_lon1 <= grid['max_lon'] and grid['min_lat'] <= chosen_lat1 <= grid['max_lat'] and np.round(grid['min_lon']) == chosen_grid1_minlon):
                        dip_val1 = grid['interp']((chosen_lat1, chosen_lon1))
                        break

                dip_val2 = None
                for grid in dip_grids:
                    if (chosen_slab_depth2 != None and grid['min_lon'] <= chosen_lon2 <= grid['max_lon'] and grid['min_lat'] <= chosen_lat2 <= grid['max_lat'] and np.round(grid['min_lon']) == chosen_grid2_minlon):
                        dip_val2 = grid['interp']((chosen_lat2, chosen_lon2))
                        break

                dip_val3 = None
                for grid in dip_grids:
                    if (chosen_slab_depth3 != None and grid['min_lon'] <= chosen_lon3 <= grid['max_lon'] and grid['min_lat'] <= chosen_lat3 <= grid['max_lat']) and np.round(grid['min_lon']) == chosen_grid3_minlon:
                        dip_val3 = grid['interp']((chosen_lat3, chosen_lon3))
                        break

                dip_val4 = None
                for grid in dip_grids:
                    if (chosen_slab_depth4 != None and grid['min_lon'] <= chosen_lon4 <= grid['max_lon'] and grid['min_lat'] <= chosen_lat4 <= grid['max_lat']) and np.round(grid['min_lon']) == chosen_grid4_minlon:
                        dip_val4 = grid['interp']((chosen_lat4, chosen_lon4))
                        break

                # now compute curvature
                shallow_K = None; deep_K = None
                dip_avg = None; dip_avg_shall = None
                # deep
                if dip_val1 != None and dip_val2 != None and dip_val1 != np.nan and dip_val2 != np.nan:
                    dist = haversine(chosen_lon1, chosen_lat1, chosen_lon2, chosen_lat2) # km
                    dip_avg = (dip_val1 + dip_val2) / 2 # degrees
                    dz      = -1.0*(chosen_slab_depth2 - chosen_slab_depth1) # km
                    ds1      = np.sqrt(dist**2 + dz**2) # km
                    ds2       = dz / np.sin(np.radians(dip_avg)) # km
                    # ds1 should equal ds2
                    K1 = (np.deg2rad(dip_val2) - np.deg2rad(dip_val1)) / ds1
                    K2 = (np.deg2rad(dip_val2) - np.deg2rad(dip_val1)) / ds2
                    deep_K = (K1 + K2) / 2
                
                # shallow
                if dip_val3 != None and dip_val4 != None and dip_val3 != np.nan and dip_val4 != np.nan:
                    dist = haversine(chosen_lon3, chosen_lat3, chosen_lon4, chosen_lat4)
                    dip_avg_shall = (dip_val3 + dip_val4) / 2
                    dz      = -1.0*(chosen_slab_depth4 - chosen_slab_depth3)
                    ds1      = np.sqrt(dist**2 + dz**2)
                    ds2       = dz / np.sin(np.radians(dip_avg_shall))
                    K1 = (np.deg2rad(dip_val4) - np.deg2rad(dip_val3)) / ds1
                    K2 = (np.deg2rad(dip_val4) - np.deg2rad(dip_val3)) / ds2
                    shallow_K = (K1 + K2) / 2

                # now find closest vc for the segment
                min_dist = 1e9; min_vc = None
                for vc in lallemand_vc:
                    dist = haversine(lon_center, lat_center, vc[0], vc[1])
                    if dist < min_dist:
                        min_dist = dist
                        min_vc = vc[2]
                if min_dist < 250:
                        chosen_vc = min_vc/10 # cm/yr
                else:
                        chosen_vc = None

                # lon1, lat1, lon2, lat2, lon_center, lat_center, lon age, lat age, age
                # segment_locs.append((coords[j][0], coords[j][1], coords[j+1][0], coords[j+1][1],lon_center, lat_center, lon_for_age, lat_out, age))

                # seg lon, seg lat, lon dip1, lat dip1, depth dip1, dip1 (270), lon dip2, lat dip2, depth dip2, dip2 (330), lon dip3, lat dip3, depth dip3, dip3 (175), lon dip4, lat dip4, depth dip4, dip4 (225)
                segment_dips.append((lon_center, lat_center, chosen_lon1, chosen_lat1, chosen_slab_depth1, dip_val1, chosen_lon2, chosen_lat2, chosen_slab_depth2, dip_val2, \
                    chosen_lon3, chosen_lat3, chosen_slab_depth3, dip_val3, chosen_lon4, chosen_lat4, chosen_slab_depth4, dip_val4))

                # all of the main quantities
                segment_data.append((lon_center, lat_center, shallow_K, deep_K, chosen_vc, age, lon_for_age, lat_out, dip_avg, dip_avg_shall))

# segment_locs = np.array(segment_locs)    
segment_dips = np.array(segment_dips)
segment_data = np.array(segment_data)
np.savetxt("data/segment_data.txt", segment_data, fmt='%s', delimiter=',')

# plot ages and segment centers
norm = mcolors.Normalize(vmin=0, vmax=200)
x_center, y_center = m1(segment_data[:,0], segment_data[:,1])
ax1.plot(x_center, y_center, 'o', color='magenta', markersize=4)
x_out, y_out = m1(segment_data[:,6], segment_data[:,7])
sc = ax1.scatter(x_out, y_out, c=segment_data[:,5], cmap='jet', norm=norm, s=30, edgecolors='k', zorder=5)
cb1 = plt.colorbar(sc, ax=ax1, label='Oceanic Age (Ma)')

# plot the dips
norm2 = mcolors.Normalize(vmin=10, vmax=80)
x_depth1, y_depth1 = m2(segment_dips[:,2], segment_dips[:,3])
cd = ax2.scatter(x_depth1, y_depth1, c=segment_dips[:,5], cmap='jet', norm=norm2, s=40, edgecolors='k', zorder=5)
x_depth2, y_depth2 = m2(segment_dips[:,6], segment_dips[:,7])
ax2.scatter(x_depth2, y_depth2, c=segment_dips[:,9], cmap='jet', norm=norm2, s=40, edgecolors='k', zorder=5)
x_depth3, y_depth3 = m2(segment_dips[:,10], segment_dips[:,11])
ax2.scatter(x_depth3, y_depth3, c=segment_dips[:,13], cmap='jet', norm=norm2, s=40, edgecolors='k', zorder=5)
x_depth4, y_depth4 = m2(segment_dips[:,14], segment_dips[:,15])
ax2.scatter(x_depth4, y_depth4, c=segment_dips[:,17], cmap='jet', norm=norm2, s=40, edgecolors='k', zorder=5)
# cb2 = plt.colorbar(cd, label='dip')

# plot the curvature points here
vminK, vmaxK = -0.0025, 0.0025
norm2 = mcolors.Normalize(vmin=vminK, vmax=vmaxK)
num_shallow_K = 0
num_deep_K = 0
num_no_K = 0
for pt in segment_data:
    lon_center, lat_center, shallow_K, deep_K = pt[0], pt[1], pt[2], pt[3]
    x, y = m2(lon_center, lat_center)
    # Decide which curvature value to use
    if deep_K != None and not np.isnan(deep_K):
        K_val = deep_K
        edge = 'red'
        num_deep_K = num_deep_K + 1
    elif shallow_K != None and not np.isnan(shallow_K):
        K_val = shallow_K
        edge = 'gray'
        num_shallow_K = num_shallow_K + 1
    else:
        num_no_K = num_no_K + 1
        continue  # skip if both are nan
    ck = ax2.scatter(x, y, s=40, c=K_val, cmap='Spectral', norm=norm2, edgecolors=edge, linewidths=1.5, zorder=10)
ck2 = plt.colorbar(ck, ax=ax2, label='K (1/km)',extend='both')
print(num_deep_K,num_shallow_K,num_no_K)

# plot vc points here
vminvc, vmaxvc = 0, 10
norm3 = mcolors.Normalize(vmin=vminvc, vmax=vmaxvc)
for j in range(len(segment_data)):
    vc = segment_data[j,4]
    lon_plot, lat_plot = segment_data[j,6], segment_data[j,7]
    x, y = m2(lon_plot, lat_plot)
    if vc != None:
        ck2 = ax2.scatter(x, y, s=40, c=vc, cmap='Spectral', norm=norm3, edgecolors='black', linewidths=1.5, zorder=10)
    else:
        ax2.scatter(x, y, s=40, c='gray', edgecolors='black', linewidths=1.5, zorder=10)
ck3 = plt.colorbar(ck2, ax=ax2, label='vc (cm/yr)',extend='both')

# # Plot ocean age points (if desired)
# plot_ocean_age(m, "data/Muller2008/age.3.6.nc")

plt.tight_layout()
plt.show()
