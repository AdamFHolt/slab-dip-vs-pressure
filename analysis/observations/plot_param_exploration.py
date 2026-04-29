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
from functions import make_strictly_ascending, compute_DP_hs, compute_DP_pl
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


slab_visc = float(4e22)           # Pa s [ref=5e22]
alpha = float(3.28e-5)               # 1/K  [ref=3e-5]
Tm = float(1333)                  # degC [ref=1300]
diffusivity = float(8.044e-7)         # m^2/s [ref=1e-6]
plate_thick = float(88e3)         # m [ref=88e3]
crust_thick = float(7e3)         # m [ref=7e3]
cooling_model = 'plate-cooling'

# conversion factors
cmyr_to_ms = 1e-2 / 3.154e7
Ma_to_s = 1e6 * 3.154e7

plotname=''.join(['plots/DP-param-exploration.png'])
plotname_pdf=''.join(['plots/DP-param-exploration.pdf'])

# ------------------------------------------------
# --- Set up the map ---
# ------------------------------------------------
fig = plt.figure(figsize=(8,6))
G = gridspec.GridSpec(2,3)

# -----------------------------------------------
# ----------- first plot ------------------
# ------------------------------------------------

# Define the independent variables:
T_vals = np.arange(1080, 1450+10, 10)
plate_thick_vals = np.arange(75e3, 135e3+2e3, 2e3)

# orginal exploration (crust thick = 7 km)
DPmean_thresh_grid = np.zeros((len(plate_thick_vals), len(T_vals)))
# Loop over all combinations of T and plate thickness and compute DPmean_thresh.
for i, plate_thick_val in enumerate(plate_thick_vals):
    for j, T_val in enumerate(T_vals):
        data = load_data_file(slab_visc, alpha, float(T_val), diffusivity, float(plate_thick_val), crust_thick, cooling_model)
        DPmean_thresh, DPmean_tot = stats_DP(data,thresh=5)
        DPmean_thresh_grid[i, j] = DPmean_thresh

# orginal exploration (crust thick = 0 km)
DPmean_thresh_grid_nocrust = np.zeros((len(plate_thick_vals), len(T_vals)))
# Loop over all combinations of T and plate thickness and compute DPmean_thresh.
for i, plate_thick_val in enumerate(plate_thick_vals):
    for j, T_val in enumerate(T_vals):
        data = load_data_file(slab_visc, alpha, float(T_val), diffusivity, float(plate_thick_val), 0.0, cooling_model)
        DPmean_thresh, DPmean_tot = stats_DP(data,thresh=5)
        DPmean_thresh_grid_nocrust[i, j] = DPmean_thresh

plate_thick_vals_km = plate_thick_vals / 1e3

# Create the subplot (ax3) for the colored grid.
ax3 = plt.subplot(G[0, 0])

contour_levels1 = np.arange(20, 40+2, 2)
contour_levels1_minor = np.arange(21, 40, 2)
contours1 = ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid, levels=contour_levels1, colors='black', linewidths=1, linestyles='solid')
labels1 = ax3.clabel(contours1, levels=contour_levels1, inline=False, fontsize=8, fmt='%d')
for txt in labels1:
    txt.set_bbox(dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid, levels=contour_levels1_minor, colors='black', linewidths=0.5, linestyles='dashed')

contour_levels2 = np.arange(20, 40+2, 2)
contour_levels2_minor = np.arange(21, 40, 2)
contours2 = ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid_nocrust, levels=contour_levels2, colors='red', linewidths=1, linestyles='solid')
labels2 = ax3.clabel(contours2, levels=contour_levels2, inline=False, fontsize=8, fmt='%d')
for txt in labels2:
    txt.set_bbox(dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid_nocrust, levels=contour_levels2_minor, colors='red', linewidths=0.5, linestyles='dashed')

# Force right-side labels on lines that exit through the bottom before reaching T=1450
def _rightmost_label(ax, cs, level, xlim, ylim, color, fontsize=8):
    lidx = np.where(np.isclose(cs.levels, level))[0]
    if len(lidx) == 0:
        return
    best_x, best_y = -np.inf, None
    for seg in cs.allsegs[lidx[0]]:
        mask = (seg[:,0] >= xlim[0]) & (seg[:,0] <= xlim[1]) & \
               (seg[:,1] >= ylim[0]) & (seg[:,1] <= ylim[1])
        if mask.any():
            pts = seg[mask]
            i = np.argmax(pts[:,0])
            if pts[i,0] > best_x:
                best_x, best_y = pts[i,0], pts[i,1]
    if best_y is not None:
        ax.text(best_x, best_y, str(int(level)), fontsize=fontsize, color=color,
                ha='right', va='bottom',
                bbox=dict(boxstyle='square,pad=0.1', fc='white', ec='none'))

_rightmost_label(ax3, contours1, 30, (1250, 1450), (80, 135), 'black')
_rightmost_label(ax3, contours1, 32, (1250, 1450), (80, 135), 'black')

# plot paramters of the map
ax3.plot(Tm, plate_thick/1e3, 'k*', markersize=9)

# Set axis labels and ticks.
ax3.set_xlabel('Basal temperature  [°C]')
ax3.set_ylabel('Max plate thickness  [km]')
ax3.set_xlim(1250, 1450)
ax3.set_ylim(80, 135)
ax3.set_xticks(np.arange(1250, 1450+50, 50))
ax3.set_yticks(np.arange(80, 131, 10))
ax3.tick_params(axis='both', labelsize=9)
ax3.labelsize = 12


# -----------------------------------------------
# ----------- second plot -----------------------
# -----------------------------------------------
k_vals = np.arange(5e-7, 1.8e-6+0.25e-7, 0.25e-7)
alpha_vals = np.arange(2.7e-5, 4.3e-5, 0.025e-5)


# orginal exploration (crust thick = 7 km)
DPmean_thresh_grid2 = np.zeros((len(k_vals), len(alpha_vals)))

# # Loop over all combinations of k and alpha
for i, k_val in enumerate(k_vals):
    for j, alpha_val in enumerate(alpha_vals):
        alpha_val = round(alpha_val, 8)
        k_val = round(k_val, 9)
        data = load_data_file(slab_visc, alpha_val, Tm, k_val, plate_thick, crust_thick, cooling_model)
        DPmean_thresh, DPmean_tot = stats_DP(data, thresh=5)
        DPmean_thresh_grid2[i, j] = DPmean_thresh


# Create the subplot (ax3) for the colored grid.
ax4 = plt.subplot(G[0, 1])

k_vals_e6 = k_vals/1e-6
alpha_vals_e5 = alpha_vals/1e-5

contour_levels3 = np.arange(20, 44+2, 2)
contour_levels3_minor = np.arange(21, 44, 2)
contours3 = ax4.contour(alpha_vals_e5, k_vals_e6, DPmean_thresh_grid2, levels=contour_levels3, colors='black', linewidths=1, linestyles='solid')
labels3 = ax4.clabel(contours3, levels=contour_levels3, inline=False, fontsize=8, fmt='%d')
for txt in labels3:
    txt.set_bbox(dict(boxstyle='square,pad=0.1', fc='white', ec='none'))
ax4.contour(alpha_vals_e5, k_vals_e6, DPmean_thresh_grid2, levels=contour_levels3_minor, colors='black', linewidths=0.5, linestyles='dashed')

_rightmost_label(ax4, contours3, 30, (2.7, 3.8), (0.7, 1.3), 'black')
_rightmost_label(ax4, contours3, 32, (2.7, 3.8), (0.7, 1.3), 'black')

# plot paramters of the map
ax4.plot(alpha/1e-5, diffusivity/1e-6, 'k*', markersize=9)

# Set axis labels and ticks.
ax4.set_xlabel(r'$\alpha$    [$10^{-5}$ K$^{-1}$]')
ax4.set_ylabel(r'$\kappa$    [$10^{-6}$ m$^{2}s^{-1}$]')
ax4.set_xlim(2.7, 3.8)
ax4.set_ylim(0.7, 1.3)
ax4.tick_params(axis='both', labelsize=9)
ax4.labelsize = 12


plt.tight_layout()
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=700)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()
