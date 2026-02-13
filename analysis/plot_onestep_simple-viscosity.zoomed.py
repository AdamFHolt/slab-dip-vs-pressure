#!/bin/python
"""
Zoomed evolution plot for ASPECT CSV outputs.

Changes vs original:
- Adds CLI arg: x_center_km (center of zoom window)
- Zoom window: 1000 km wide, 0–1000 km depth
- Crops interpolated grids before plotting (faster + cleaner)
- Clamps zoom window to domain bounds
- Slightly smaller figure size (still high DPI output)
- More robust font handling (falls back if Myriad not found)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os
from functions import create_grid, get_slab_midplane

# --------------------------- style / fonts ---------------------------
try:
    import matplotlib.font_manager as fm
    font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
    if os.path.exists(font_path):
        myriad_pro = fm.FontProperties(fname=font_path)
        mpl.rcParams['font.family'] = 'Myriad Pro'
    else:
        myriad_pro = None
except Exception:
    myriad_pro = None

mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['axes.labelpad'] = 1.25
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['xtick.major.pad'] = 2
mpl.rcParams['ytick.major.pad'] = 2
mpl.rcParams['xtick.major.size'] = 2.5
mpl.rcParams['ytick.major.size'] = 2.5
mpl.rcParams['xtick.minor.size'] = 1.25
mpl.rcParams['ytick.minor.size'] = 1.25

# --------------------------- CLI ---------------------------
if len(sys.argv) < 4:
    print("Usage: plot_zoom.py MODEL_NAME TIME_STEP X_CENTER_KM")
    print("Example: plot_zoom.py lhr_caseA 120 3000")
    sys.exit(1)

model_name = str(sys.argv[1])
time = int(sys.argv[2])
x_center_km = float(sys.argv[3])  # NEW: center of zoom window

# --------------------------- analysis depths ---------------------------
analysis_depth1 = 250e3
analysis_depth2 = 300e3
analysis_depth3 = 350e3

# --------------------------- model properties ---------------------------
xmax = 5800.e3
ymax = 1450.e3

# --------------------------- ASPECT output ---------------------------
csvs_loc = 'csv_outputs/'
models_loc = 'raw_outputs/'
stats_file = ''.join([models_loc, str(model_name), '/statistics'])
model_output_dt = 50      # output dt as set in ASPECT .prm file
num_header_lines = 16     # num header lines in stats_file

out_dir = 'plots/evolution/simple/zoomed/'
os.makedirs(out_dir, exist_ok=True)

plotname = ''.join([out_dir, str(model_name), '.', str(time), '.png'])
plotname_pdf = ''.join([out_dir, str(model_name), '.', str(time), '.pdf'])

# --------------------------- CSV columns ---------------------------
c_llith_col = 25
P_col = 29
x_col = 30
y_col = 31
visc_col = 26
vx_col = 0
vz_col = 1

# --------------------------- grids ---------------------------
xmin_plot = 0
ymin_plot = 0
grid_res = 2.0e3

X_low, Y_low = create_grid(xmin_plot, xmax, ymin_plot, ymax, grid_res)

# extremely low-res grid for velocities
X_vels, Y_vels = create_grid(xmin_plot, xmax, ymin_plot, ymax, 10.e3)

# higher-res grid for detailed calculations (slab contour extraction)
xmin_plot2 = 2500.e3
xmax_plot2 = 3750.e3
ymin_plot2 = ymax - 600.e3
grid_res2 = 1.0e3
X_low2, Y_low2 = create_grid(xmin_plot2, xmax_plot2, ymin_plot2, ymax, grid_res2)

# --------------------------- zoom window (km) ---------------------------
zoom_w_km = 1200.0
zoom_d_km = 1000.0

# Convert to meters for clamping; but we use km for plotting limits.
x0_km = x_center_km - zoom_w_km / 2.0
x1_km = x_center_km + zoom_w_km / 2.0
z0_km = 0.0
z1_km = zoom_d_km

# Clamp to domain in km
domain_x0_km = 0.0
domain_x1_km = xmax / 1.e3
x0_km = max(domain_x0_km, x0_km)
x1_km = min(domain_x1_km, x1_km)

# If window got clipped at an edge, re-expand to keep width if possible
if (x1_km - x0_km) < zoom_w_km:
    if x0_km == domain_x0_km:
        x1_km = min(domain_x1_km, x0_km + zoom_w_km)
    elif x1_km == domain_x1_km:
        x0_km = max(domain_x0_km, x1_km - zoom_w_km)

# --------------------------- load CSV + dimensional time ---------------------------
csv_filename = ''.join([csvs_loc, model_name, '/full.', str(time), '.csv'])

stats_line_num = num_header_lines + (time * model_output_dt)
with open(stats_file, 'r') as f:
    line = f.readlines()[stats_line_num]
time_dim = float(line.split()[1]) / 1.e6  # Myr
print("%.0f: t = %.1f Myr" % (time, time_dim))

model_data = np.loadtxt(csv_filename, delimiter=',', skiprows=1)

# --------------------------- interpolate fields ---------------------------
visc = griddata(
    (model_data[:, x_col], model_data[:, y_col]),
    model_data[:, visc_col],
    (X_low, Y_low),
    method='nearest'
)

llith = griddata(
    (model_data[:, x_col], model_data[:, y_col]),
    model_data[:, c_llith_col],
    (X_low2, Y_low2),
    method='linear'
)

vx = griddata(
    (model_data[:, x_col], model_data[:, y_col]),
    model_data[:, vx_col],
    (X_vels, Y_vels),
    method='nearest'
)

vz = griddata(
    (model_data[:, x_col], model_data[:, y_col]),
    model_data[:, vz_col],
    (X_vels, Y_vels),
    method='nearest'
)

vmag = np.sqrt(vx**2 + vz**2)

# --------------------------- slab midplane points (km) ---------------------------
# comp_contour_val = 0.5
# llith_cont = plt.contour(X_low2/1.e3, (ymax - Y_low2)/1.e3, llith, levels=[comp_contour_val])
# llith_points_tmp = llith_cont.allsegs[0][0]

# cutoff_shall = 110.0
# cutoff_deep = 575.0
# llith_points = get_slab_midplane(llith_points_tmp, cutoff_shall, cutoff_deep)

# def get_midpoint_at_depth(llith_pts, depth_km):
#     misfit = 1e9
#     xcent = np.nan
#     zcent = np.nan
#     for j in range(len(llith_pts)):
#         misfit_tmp = np.abs(llith_pts[j, 1] - depth_km)
#         if misfit_tmp < misfit:
#             misfit = misfit_tmp
#             xcent = llith_pts[j, 0]
#             zcent = llith_pts[j, 1]
#     return xcent, zcent

# xcent1, zcent1 = get_midpoint_at_depth(llith_points, analysis_depth1/1.e3)
# xcent2, zcent2 = get_midpoint_at_depth(llith_points, analysis_depth2/1.e3)
# xcent3, zcent3 = get_midpoint_at_depth(llith_points, analysis_depth3/1.e3)

# --------------------------- crop arrays to zoom window ---------------------------
# X_low, Y_low are 2D; extract 1D axes from a row/col.
x_low_km_1d = X_low[0, :] / 1.e3
z_low_km_1d = (ymax - Y_low[:, 0]) / 1.e3  # depth km, 0 at surface

ix_low = np.where((x_low_km_1d >= x0_km) & (x_low_km_1d <= x1_km))[0]
iz_low = np.where((z_low_km_1d >= z0_km) & (z_low_km_1d <= z1_km))[0]

if len(ix_low) < 2 or len(iz_low) < 2:
    raise RuntimeError("Zoom window too small / out of bounds for X_low/Y_low grid.")

X_low_z = X_low[np.ix_(iz_low, ix_low)]
Y_low_z = Y_low[np.ix_(iz_low, ix_low)]
visc_z = visc[np.ix_(iz_low, ix_low)]

# velocity grid crop
x_vel_km_1d = X_vels[0, :] / 1.e3
z_vel_km_1d = (ymax - Y_vels[:, 0]) / 1.e3

ix_vel = np.where((x_vel_km_1d >= x0_km) & (x_vel_km_1d <= x1_km))[0]
iz_vel = np.where((z_vel_km_1d >= z0_km) & (z_vel_km_1d <= z1_km))[0]

if len(ix_vel) < 2 or len(iz_vel) < 2:
    raise RuntimeError("Zoom window too small / out of bounds for X_vels/Y_vels grid.")

X_vels_z = X_vels[np.ix_(iz_vel, ix_vel)]
Y_vels_z = Y_vels[np.ix_(iz_vel, ix_vel)]
vx_z = vx[np.ix_(iz_vel, ix_vel)]
vz_z = vz[np.ix_(iz_vel, ix_vel)]
vmag_z = vmag[np.ix_(iz_vel, ix_vel)]

# --------------------------- plotting ---------------------------
fig = plt.figure(figsize=(3.4, 3.4))  # smaller plot (square-ish)
gs = GridSpec(1, 1)

ax1 = fig.add_subplot(gs[0, 0])

# viscosity field
visc_plot = ax1.contourf(
    X_low_z/1.e3,
    (ymax - Y_low_z)/1.e3,
    np.log10(visc_z),
    cmap=cm.get_cmap('hot_r'),
    levels=np.linspace(20, 25, 501),
    extend='max'
)

# flow vectors
ax1.streamplot(
    np.flipud(X_vels_z/1.e3),
    np.flipud((ymax - Y_vels_z)/1.e3),
    np.flipud(vx_z),
    np.flipud(-1.0 * vz_z),
    color='dimgray',
    linewidth=vmag_z * 18,
    arrowsize=0.3,
    density=1.2
)

# # colorbar (tight to the right)
# cbar = plt.colorbar(
#     visc_plot,
#     cax=fig.add_axes([0.88, 0.18, 0.03, 0.64]),
#     ticks=[20, 21, 22, 23, 24, 25],
#     ticklocation='right'
# )
# cbar.ax.tick_params(axis='y', labelsize=5.5, pad=1,
#                     left=False, labelleft=False, right=True, labelright=True)

# axis properties: zoom box
ax1.set_xlim([x0_km, x1_km])
ax1.set_ylim([z1_km, z0_km])  # invert so depth increases downward
ax1.set_aspect('equal', adjustable='box')
ax1.tick_params(direction='out', length=2, labelsize=6)

# time label
ax1.annotate(
    f"{time_dim:.1f} Myr",
    xy=(0.98, 1.04),
    xycoords='axes fraction',
    va='center',
    ha='right',
    fontsize=10,
    color='k'
)

# # mark analysis point (only if it's inside zoom box)
# if (xcent2 >= x0_km) and (xcent2 <= x1_km) and (zcent2 >= z0_km) and (zcent2 <= z1_km):
#     ax1.scatter(xcent2, zcent2, s=6, color='black', zorder=5)

# optional: label axes (comment out if you want super-clean)
ax1.set_xlabel("x (km)")
ax1.set_ylabel("depth (km)")

# save
plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
plt.savefig(plotname_pdf, bbox_inches='tight', format='pdf')
plt.close(fig)
