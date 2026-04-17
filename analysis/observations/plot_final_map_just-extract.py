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


# constants
slab_visc = float(sys.argv[1])           # Pa s [ref=5e22]
alpha = float(sys.argv[2])               # 1/K  [ref=3e-5]
Tm = float(sys.argv[3])                  # degC [ref=1300]
diffusivity = float(sys.argv[4])         # m^2/s [ref=1e-6]
plate_thick = float(sys.argv[5])         # m [ref=88e3]
crust_thick = float(sys.argv[6])         # m [ref=7e3]
hs_or_pl = 2                             # = 1 for hs and 2 for pl

# conversion factors
cmyr_to_ms = 1e-2 / 3.154e7
Ma_to_s = 1e6 * 3.154e7

# plot/file names
if hs_or_pl == 1:
    cooling_model = 'hs-cooling'
else:
    cooling_model = 'plate-cooling'

textname=''.join(['text_files/new/maps.slab',str(slab_visc),'.alpha',str(alpha),'.T',str(Tm),'.k',str(diffusivity),'.platethick',str(plate_thick),'.crustthick',str(crust_thick),'.',cooling_model,'.txt'])

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
DPmin, DPmax = 0, 60
vminK, vmaxK = 0, 20
norm1 = mcolors.Normalize(vmin=DPmin, vmax=DPmax)
norm2 = mcolors.Normalize(vmin=vminK, vmax=vmaxK)
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


    if not np.isnan(dip_shall) and age < 250 and not np.isnan(vc) and not np.isnan(K):

        # compute/plot ηHKvc/L_eff -------------------
        H_eff = compute_H_eff(age, Tm, k=diffusivity, plate_thick=plate_thick)  # m
        stress_scaling = K * (vc * cmyr_to_ms) * slab_visc * H_eff / 1.497e6 * 1e-6  # MPa

        # compute/plot DP ------------------------
        if not np.isnan(dip_deep):
            if hs_or_pl == 1:
                DP = compute_DP_hs(age, dip_deep, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick)
            else:
                DP = compute_DP_pl(age, dip_deep, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
        else:
            if hs_or_pl == 1:
                DP = compute_DP_hs(age, dip_shall, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick)
            else:
                DP = compute_DP_pl(age, dip_shall, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
                
        # store DP and stress scaling
        DP_array.append((DP, np.abs(stress_scaling)))

DP_array = np.array(DP_array)
np.savetxt(textname, DP_array, delimiter=',', fmt='%f')
