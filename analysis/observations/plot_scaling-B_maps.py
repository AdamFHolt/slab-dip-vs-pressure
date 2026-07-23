#!/usr/bin/env python
"""
Supplementary figure: breakdown of the nondimensional scaling parameter
Lambda = [(eta*H*|dtheta/ds|*v_c)/L_eff] / B into its two ingredients:
  (a) the dimensional shear-stress scaling (eta*H*|dtheta/ds|*v_c)/L_eff [MPa]
  (b) the depth-integrated slab buoyancy B = Drho*g*H [MPa] (no cos(theta))
Both use the reference parameter set (Richards et al., 2018 "Pk" model with
T fixed at 1333 C) and reference slab viscosity 4e22 Pa s.
"""
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as mcolors
from matplotlib.patches import PathPatch
import numpy as np
from netCDF4 import Dataset
import glob
from functions import read_pb2002_boundaries
from functions import make_strictly_ascending, compute_H_eff, compute_B_pl
from cmcrameri import cm as cmc
import matplotlib.font_manager as fm
font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
myriad_pro = fm.FontProperties(fname=font_path)


mpl.rcParams['font.family'] = 'Myriad Pro'
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

# reference parameters (Richards et al. 2018, simple plate model, T=1333C)
slab_visc = 4e22        # Pa s
Tm = 1333.0             # degC
diffusivity = 8.044e-7  # m^2/s
alpha = 3.28e-5         # 1/K
plate_thick = 88e3      # m
crust_thick = 7e3       # m

cmyr_to_ms = 1e-2 / 3.154e7

plotname = 'plots/maps_scaling-B.png'
plotname_pdf = 'plots/maps_scaling-B.pdf'

# ------------------------------------------------
# --- Set up the maps ---
# ------------------------------------------------
fig = plt.figure(figsize=(8,6))
G = gridspec.GridSpec(2,3)

ax1 = plt.subplot(G[0, 0:2])
m1 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax1)
m1.drawcoastlines(linewidth=0)
limb1 = m1.drawmapboundary(fill_color='white')
for _poly in m1.fillcontinents(color='silver', lake_color='silver'):
    _poly.set_clip_path(limb1.get_path(), limb1.get_transform())
ax1.add_patch(PathPatch(limb1.get_path(), transform=limb1.get_transform(),
                          facecolor='none', edgecolor='black', linewidth=1, zorder=50))

ax2 = plt.subplot(G[1, 0:2])
m2 = Basemap(projection='robin', lon_0=-180, resolution='l',
            llcrnrlat=-80, urcrnrlat=80, ax=ax2)
m2.drawcoastlines(linewidth=0)
limb2 = m2.drawmapboundary(fill_color='white')
for _poly in m2.fillcontinents(color='silver', lake_color='silver'):
    _poly.set_clip_path(limb2.get_path(), limb2.get_transform())
ax2.add_patch(PathPatch(limb2.get_path(), transform=limb2.get_transform(),
                          facecolor='none', edgecolor='black', linewidth=1, zorder=50))

#------------------------------------------------
# ------ plot slab contours ---------------------
#------------------------------------------------
slab_files = glob.glob("data/Slab2/*_dep_*grd")
for grd_file in slab_files:
    nc = Dataset(grd_file, 'r')
    lon_arr = nc.variables['x'][:]
    lat_arr = nc.variables['y'][:]
    slab_data = nc.variables['z'][:]
    lat_arr = make_strictly_ascending(lat_arr)
    lon_arr = make_strictly_ascending(lon_arr)
    nc.close()

    lon2d, lat2d = np.meshgrid(lon_arr, lat_arr)
    levels = np.arange(-600, 1, 100)
    x2d, y2d = m1(lon2d, lat2d)
    ax1.contour(x2d, y2d, slab_data, levels=levels, colors='red', linestyles='solid', linewidths=0.5, zorder=10)
    x2d2, y2d2 = m2(lon2d, lat2d)
    ax2.contour(x2d2, y2d2, slab_data, levels=levels, colors='red', linestyles='solid', linewidths=0.5, zorder=10)

# ------------------------------------------------
# -------- Read in the segment data --------------
# ------------------------------------------------
data = np.genfromtxt("data/segment_data.txt", delimiter=',', dtype=str)
data[data == 'None'] = np.nan #  lon_center, lat_center, shallow_K, deep_K, chosen_vc, age, lon_for_age, lat_out, dip, dip_shall
segment_data = np.array(data, dtype=float)

# ------------------------------------------------
# ---- compute and plot scaling and B ------------
# ------------------------------------------------
norm1 = mcolors.Normalize(vmin=0, vmax=15)   # stress scaling [MPa]
norm2 = mcolors.Normalize(vmin=30, vmax=80)  # B [MPa]
for i in range(len(segment_data)):

    lon_center, lat_center = segment_data[i,0], segment_data[i,1]
    age = segment_data[i,5]               # Ma
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

        x1, y1 = m1(lon_center, lat_center)

        # Lambda numerator: shear-stress scaling [MPa]
        H_eff = compute_H_eff(age, Tm, k=diffusivity, plate_thick=plate_thick)  # m
        stress_scaling = K * (vc * cmyr_to_ms) * slab_visc * H_eff / 1.497e6 * 1e-6  # MPa
        cs = ax1.scatter(x1, y1, s=23, c=np.abs(stress_scaling), cmap=cmc.navia_r, norm=norm1, edgecolors='black', linewidths=0.4, zorder=10)

        # Lambda denominator: net buoyancy B = Drho*g*H (no cos(dip)) [MPa]
        B_seg = compute_B_pl(age, Tm, k=diffusivity, rho0=3330., alpha=alpha, crust_density=3450, crust_thick=crust_thick, plate_thick=plate_thick)
        # navia (un-reversed): low B = dark, matching its effect on Lambda
        cb = ax2.scatter(x1, y1, s=23, c=B_seg, cmap=cmc.navia, norm=norm2, edgecolors='black', linewidths=0.4, zorder=10)

cs_bar = plt.colorbar(cs, ax=ax1, extend='max', shrink=0.5, pad=0.05)
cs_bar.set_ticks([0, 5, 10, 15])
cs_bar.ax.set_title(r'$(\eta H |d\theta/ds|\, v_c)/L_\mathrm{eff}$  [MPa]', fontsize=10, pad=8)
cb_bar = plt.colorbar(cb, ax=ax2, extend='min', shrink=0.5, pad=0.05)
cb_bar.set_ticks([30, 40, 50, 60, 70, 80])
cb_bar.ax.set_title(r'$B$  [MPa]', fontsize=10, pad=8)

# ------------------------------------------------
# --- Plot plate boundaries ---
# ------------------------------------------------
boundaries = read_pb2002_boundaries("data/PB2002/PB2002_boundaries.dig.txt")
for seg in boundaries:
    if len(seg['coords']) < 2:
        continue
    lons, lats = zip(*seg['coords'])
    x, y = m1(lons, lats)
    if seg['symbol'] in ["/", "\\"]:
        color = 'magenta'
        lw = 1.25
    else:
        color = 'black'
        lw = 0.5
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
