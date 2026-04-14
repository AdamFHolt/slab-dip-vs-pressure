#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from netCDF4 import Dataset
import glob
import math
import scipy.special
import scipy.integrate


def make_strictly_ascending(arr, eps=1e-8):
    """
    Ensure the array is strictly ascending by adding a small epsilon where needed.
    """
    arr = np.array(arr, copy=True)
    for i in range(1, len(arr)):
        if arr[i] <= arr[i - 1]:
            arr[i] = arr[i - 1] + eps
    return arr


def read_pb2002_boundaries(filename):
    """
    Read plate boundaries from a PB2002 text file.
    Each segment begins with a header line whose first token is 5 characters:
      - characters 0-1: left plate ID
      - character 2: boundary symbol ("/", "\" for subduction; "-" for non-subduction)
      - characters 3-4: right plate ID
    Subsequent lines contain coordinate pairs (lon, lat) separated by a comma.
    A line starting with "***" indicates the end of a segment.
    Returns a list of segments. Each segment is a dict with keys:
      'left_plate', 'symbol', 'right_plate', and 'coords' (a list of (lon, lat) tuples).
    """
    segments = []
    current_segment = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Check for end-of-segment marker.
            if line.startswith("***"):
                current_segment = None
                continue
            # Split line into tokens.
            tokens = line.split()
            # If the first token is exactly 5 characters and matches expected pattern, assume header.
            if tokens and len(tokens[0]) == 5 and tokens[0][0:2].isalpha() and tokens[0][2] in "/-\\" and tokens[0][3:5].isalpha():
                header = tokens[0]
                left_plate = header[0:2]
                symbol = header[2]
                right_plate = header[3:5]
                current_segment = {
                    'left_plate': left_plate,
                    'symbol': symbol,
                    'right_plate': right_plate,
                    'coords': []
                }
                segments.append(current_segment)
            else:
                # Otherwise, assume it's a coordinate line.
                try:
                    # Expecting a format like: "+7.77235E+00,-5.43960E+01"
                    parts = line.split(',')
                    if len(parts) == 2:
                        lon = float(parts[0])
                        lat = float(parts[1])
                        if current_segment is not None:
                            current_segment['coords'].append((lon, lat))
                except Exception as e:
                    print("Error parsing coordinate line:", line, e)
    return segments


def plot_ocean_age(m,filename):
    # ---------------------------
    ocean_age_file = filename
    nc = Dataset(ocean_age_file, 'r')
    lon_arr = nc.variables['x'][:]  # Longitude
    lat_arr = nc.variables['y'][:]  # Latitude
    age_data = nc.variables['z'][:]  # Oceanic age data

    try:
        age_data, lon_arr = shiftgrid(0., age_data, lon_arr, start=False)
    except Exception as e:
        lon_arr = np.where(lon_arr > 180, lon_arr - 360, lon_arr)

    # Create 2D meshgrid for longitude & latitude
    lon2d, lat2d = np.meshgrid(lon_arr, lat_arr)

    # Project to map coordinates
    x2d, y2d = m(lon2d, lat2d)


def haversine(lon1, lat1, lon2, lat2):
    """
    Compute the great-circle distance between two points on Earth in km.
    """
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    return 6371 * c

def calculate_bearing(lon1, lat1, lon2, lat2):
    """
    Calculate the initial bearing (forward azimuth) from point 1 to point 2.
    """
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    initial_bearing = math.degrees(math.atan2(x, y))
    return (initial_bearing + 360) % 360

def destination_point(lon, lat, bearing, distance_km):
    """
    Given a start point (lon, lat), a bearing (in degrees), and a distance (km),
    return the destination point using spherical trigonometry.
    """
    R = 6371.0  # Earth's radius in km
    bearing_rad = math.radians(bearing)
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)
    
    lat2 = math.asin(math.sin(lat_rad) * math.cos(distance_km/R) +
                     math.cos(lat_rad) * math.sin(distance_km/R) * math.cos(bearing_rad))
    lon2 = lon_rad + math.atan2(math.sin(bearing_rad) * math.sin(distance_km/R) * math.cos(lat_rad),
                                math.cos(distance_km/R) - math.sin(lat_rad) * math.sin(lat2))
    return math.degrees(lon2), math.degrees(lat2)



def load_data_file(slab_visc,alpha,Tm,diffusivity,plate_thick,crust_thick,cooling_model):

    name=''.join(['text_files/new/maps.slab',str(slab_visc),'.alpha',str(alpha),'.T',str(Tm),'.k',str(diffusivity),'.platethick',str(plate_thick),'.crustthick',str(crust_thick),'.',cooling_model,'.txt'])

    array = np.genfromtxt(name, delimiter=',', dtype=float)

    return array


def stats_data_file(array,thresh):

    sum_tot = 0;  n_tot = 0; n_perc = 0

    for i in range(len(array)):

        if not np.isnan(array[i,1]):
            sum_tot += array[i,1]
            n_tot += 1
        
        if np.abs(array[i,1]) < thresh:
            n_perc += 1
    
    perc = (n_perc/n_tot) *100
    mean = sum_tot/n_tot
    std  = np.std(array[:,1])

    return mean, std, perc


def stats_DP(array,thresh):

    DPtot_all = 0; DPtot_thresh = 0; 
    n_all = 0; n_thresh = 0

    for i in range(len(array)):

        DP = array[i,0]
        scaling = array[i,1]

        DPtot_all += DP
        n_all += 1

        if np.abs(scaling) < thresh:
            DPtot_thresh += DP
            n_thresh += 1
    
    DPmean_all = DPtot_all/n_all
    DPmean_thresh = DPtot_thresh/n_thresh

    return DPmean_thresh, DPmean_all


def compute_DP_hs(age, dip, Tm, k, rho0, alpha, crust_density, crust_thick):

    # DT=1300.		# K
    # k=1.e-6		# m2/s
    # rho0=3300.	# kg/m3
    # alpha=3.e-5 	# 1/K
    # crust? (if no, set crust_thick = 0)
    # crust_density  = 3450	# rho_eclogite [kg/m3], C-T Lee, W-P Chen, 2007, EPSL
    # crust_thick = 0	# m

    g=9.81		# m/s2

    # set up depths, compute T and density, integrate density
    z=np.arange(0,400.0,0.25)*1e3 # m
    T_erf= Tm - Tm * scipy.special.erfc(z/(2*np.sqrt(k*age*1e6*365*24*60*60)))

    B_erf= (Tm - T_erf) * rho0 * alpha # rho - rho0, kg/m3
    B_erf_int=scipy.integrate.simpson(y=B_erf, x=z) + (crust_thick * (crust_density-rho0)) # kg/m2
    DP = np.cos(np.deg2rad(dip))*B_erf_int*g*1e-6 # MPa

    return DP


def compute_H_eff(age, Tm, k, plate_thick, frac=0.9):
    """
    Return the depth (m) at which the plate cooling temperature reaches frac*Tm.

    Uses the same plate cooling series as compute_DP_pl.  T increases from 0
    at the surface (z=0) to Tm at z=plate_thick; H_eff is the depth at which
    the plate has warmed to frac*Tm (i.e., where thermal thickness ends).

    Parameters
    ----------
    age        : float, plate age in Ma
    Tm         : float, mantle temperature in K (or degC — only the ratio matters)
    k          : float, thermal diffusivity in m^2/s
    plate_thick: float, asymptotic plate thickness in m
    frac       : float, fraction of Tm that defines the base of the plate (default 0.9)

    Returns
    -------
    H_eff : float, depth in m where T_plate = frac*Tm; returns plate_thick if
            the temperature never reaches frac*Tm within the plate.
    """
    age_s = age * 1e6 * 365 * 24 * 60 * 60  # Ma -> s
    z = np.arange(0, plate_thick / 1e3, 0.1) * 1e3  # m

    T_term1 = Tm * (z / plate_thick)
    T_term2 = np.zeros_like(z)
    for n in range(1, 10):
        T_term2 += ((2 * Tm) / (n * np.pi)) * np.sin(n * np.pi * z / plate_thick) * \
                   np.exp((-1. * n**2 * np.pi**2 * k * age_s) / plate_thick**2)
    T_plate = T_term1 + T_term2

    T_target = frac * Tm
    # Find first depth where T_plate crosses T_target
    idx = np.where(T_plate >= T_target)[0]
    if len(idx) == 0:
        return plate_thick
    # Linear interpolation between the bracketing points
    i = idx[0]
    if i == 0:
        return z[0]
    z0, z1 = z[i - 1], z[i]
    T0, T1 = T_plate[i - 1], T_plate[i]
    H_eff = z0 + (T_target - T0) / (T1 - T0) * (z1 - z0)
    return H_eff


def compute_DP_pl(age, dip, Tm, k, rho0, alpha, crust_density, crust_thick, plate_thick):

    g=9.81		# m/s2
    age = age * 1e6 * 365 * 24 * 60 * 60 # seconds

    # set up depths, compute T and density, integrate density
    z=np.arange(0,plate_thick/1e3,0.1)*1e3 # m

    T_term1 = Tm*(z/plate_thick)
    T_term2 = 0
    for n in range(1,10):
        T_term2 += ((2*Tm)/(n*np.pi))*np.sin(n*np.pi*z/plate_thick)*np.exp((-1.* n**2 * (np.pi)**2 * k * age)/(plate_thick**2))
    T_plate = T_term1 + T_term2

    B_plate= (Tm - T_plate) * rho0 * alpha # rho - rho0, kg/m3

    B_plate_int=scipy.integrate.simpson(y=B_plate, x=z) + (crust_thick * (crust_density-rho0)) # kg/m2
    DP = np.cos(np.deg2rad(dip))*B_plate_int*g*1e-6 # MPa


    return DP
