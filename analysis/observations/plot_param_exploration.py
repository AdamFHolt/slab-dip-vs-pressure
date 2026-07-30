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
from functions import qualifying_mask, DPmean_fixed
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
lambda_thresh = 0.1               # threshold for nondimensional scaling parameter Lambda

# conversion factors
cmyr_to_ms = 1e-2 / 3.154e7
Ma_to_s = 1e6 * 3.154e7

plotname=''.join(['plots/DP-param-exploration.png'])
plotname_pdf=''.join(['plots/DP-param-exploration.pdf'])

def _label_points(cs, xlim, ylim, avoid, placed, pad=0.07):
    """One label point per contour level: the point on the level's longest
    visible run that sits farthest from the reference star and from the labels
    already placed, so the crust and no-crust families do not collide."""
    xr, yr = xlim[1] - xlim[0], ylim[1] - ylim[0]
    pts = []
    for segs in cs.allsegs:
        best = None
        for seg in segs:
            vis = seg[(seg[:, 0] >= xlim[0] + pad * xr) & (seg[:, 0] <= xlim[1] - pad * xr) &
                      (seg[:, 1] >= ylim[0] + pad * yr) & (seg[:, 1] <= ylim[1] - pad * yr)]
            if len(vis) > 1 and (best is None or len(vis) > len(best)):
                best = vis
        if best is None:
            continue
        others = placed + [avoid]
        clear = lambda p: min(((p[0] - q[0]) / xr)**2 + ((p[1] - q[1]) / yr)**2 for q in others)
        p = max(best[:: max(1, len(best) // 40)], key=clear)
        placed.append((p[0], p[1]))
        pts.append((p[0], p[1]))
    return pts


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

# The segment population is frozen at the reference parameters (the star), so the
# sweep averages DP over the same segments everywhere and the field is continuous.
# Re-selecting on Lambda at every parameter combination instead makes the mean step
# as segments cross the threshold, which shows up as kinked contours.  The same
# reference population is used for both the crust and no-crust families, so the
# difference between them is the eclogite crust term alone.
ref_data = load_data_file(slab_visc, alpha, Tm, diffusivity, plate_thick, crust_thick, cooling_model)
ref_mask = qualifying_mask(ref_data, lambda_thresh)
print('reference population: %d of %d segments, mean DP = %.2f MPa'
      % (ref_mask.sum(), len(ref_data), DPmean_fixed(ref_data, ref_mask)))

# orginal exploration (crust thick = 7 km)
DPmean_thresh_grid = np.zeros((len(plate_thick_vals), len(T_vals)))
# Loop over all combinations of T and plate thickness and average over ref_mask.
for i, plate_thick_val in enumerate(plate_thick_vals):
    for j, T_val in enumerate(T_vals):
        data = load_data_file(slab_visc, alpha, float(T_val), diffusivity, float(plate_thick_val), crust_thick, cooling_model)
        DPmean_thresh_grid[i, j] = DPmean_fixed(data, ref_mask)

# orginal exploration (crust thick = 0 km)
DPmean_thresh_grid_nocrust = np.zeros((len(plate_thick_vals), len(T_vals)))
# Loop over all combinations of T and plate thickness and average over ref_mask.
for i, plate_thick_val in enumerate(plate_thick_vals):
    for j, T_val in enumerate(T_vals):
        data = load_data_file(slab_visc, alpha, float(T_val), diffusivity, float(plate_thick_val), 0.0, cooling_model)
        DPmean_thresh_grid_nocrust[i, j] = DPmean_fixed(data, ref_mask)

plate_thick_vals_km = plate_thick_vals / 1e3

# Create the subplot (ax3) for the colored grid.
ax3 = plt.subplot(G[0, 0])

contour_levels1 = np.arange(20, 40+2, 2)
contour_levels1_minor = np.arange(21, 40, 2)
contours1 = ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid, levels=contour_levels1, colors='black', linewidths=1, linestyles='solid')
ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid, levels=contour_levels1_minor, colors='black', linewidths=0.5, linestyles='dashed')

contour_levels2 = np.arange(20, 40+2, 2)
contour_levels2_minor = np.arange(21, 40, 2)
contours2 = ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid_nocrust, levels=contour_levels2, colors='red', linewidths=1, linestyles='solid')
ax3.contour(T_vals, plate_thick_vals_km, DPmean_thresh_grid_nocrust, levels=contour_levels2_minor, colors='red', linewidths=0.5, linestyles='dashed')

# label both families once the two sets of lines exist, keeping them apart
xlim3, ylim3 = (1250, 1450), (80, 135)
placed3 = []
ax3.clabel(contours1, inline=True, inline_spacing=2, fontsize=8, fmt='%d',
           manual=_label_points(contours1, xlim3, ylim3, (Tm, plate_thick/1e3), placed3))
ax3.clabel(contours2, inline=True, inline_spacing=2, fontsize=8, fmt='%d',
           manual=_label_points(contours2, xlim3, ylim3, (Tm, plate_thick/1e3), placed3))

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
ax3.set_box_aspect(1)
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
        DPmean_thresh_grid2[i, j] = DPmean_fixed(data, ref_mask)


# Create the subplot (ax3) for the colored grid.
ax4 = plt.subplot(G[0, 1])

k_vals_e6 = k_vals/1e-6
alpha_vals_e5 = alpha_vals/1e-5

contour_levels3 = np.arange(20, 44+2, 2)
contour_levels3_minor = np.arange(21, 44, 2)
contours3 = ax4.contour(alpha_vals_e5, k_vals_e6, DPmean_thresh_grid2, levels=contour_levels3, colors='black', linewidths=1, linestyles='solid')
ax4.contour(alpha_vals_e5, k_vals_e6, DPmean_thresh_grid2, levels=contour_levels3_minor, colors='black', linewidths=0.5, linestyles='dashed')

xlim4, ylim4 = (2.7, 3.8), (0.7, 1.3)
ax4.clabel(contours3, inline=True, inline_spacing=2, fontsize=8, fmt='%d',
           manual=_label_points(contours3, xlim4, ylim4, (alpha/1e-5, diffusivity/1e-6), []))


# plot paramters of the map
ax4.plot(alpha/1e-5, diffusivity/1e-6, 'k*', markersize=9)

# Set axis labels and ticks.
ax4.set_xlabel(r'$\alpha$    [$10^{-5}$ K$^{-1}$]')
ax4.set_ylabel(r'$\kappa$    [$10^{-6}$ m$^{2}s^{-1}$]')
ax4.set_xlim(2.7, 3.8)
ax4.set_ylim(0.7, 1.3)
ax4.tick_params(axis='both', labelsize=9)
ax4.set_box_aspect(1)
ax4.labelsize = 12


plt.tight_layout()
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=700)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.clf()
