#!/bin/python
"""
Variable plate thickness check (reviewer 2, comment 2): does the dQ/ds scaling
eta*H*K*v_c/L_eff, calibrated on H = 100 km, hold at H = 70 and 130 km?

eta' = 250, 500, 1000 crossed with H = 70, 100, 130 km, for free, fixed-SP and
fixed-OP plates (27 models), plus the B-matched free run (eta' = 500, H = 70,
drho = 71.43 so B = 49 MPa, as at H = 100).  One row of panels per plate
configuration; models not yet extracted are skipped with a note.
Reads text_files/withT_varH/ (built by make_withT.varH.py), whose H = 100 files
come from extracted_coverage/ and the others from extracted_varH/, so all nine
are from the same version of extract_properties.py.  Timesteps whose midplane
does not resolve the analysis depth (coverage flag, column 38) are dropped, as
in withT_covered/.

Panels follow Figs 5-6 of the manuscript, all in MPa:
  (a) DP against B_slab (the sensor test, Fig. 5)
  (b) dQ/ds against the scaling (the calibration, Fig. 6b)
  (c) B_slab - DP against the scaling (the application, Fig. 6c)

Writes plots/DP-comparisons/compilations/forces-vs-scaling.varH.z<d>...png/pdf
and prints, per model and per H, the correlation and the implied L_eff
(median of eta*H*K*v_c / (dQ/ds) over timesteps where both are positive, the
same sign filter as the L_eff = 1497 km calibration), per plate configuration
and pooled.

Optional fifth and sixth arguments choose the scaling's H exponent n and its
length scale: n = 1 is the paper's eta*H*K*v_c/L_eff with L_eff = 1497 km;
n = 3 is the moment-balance form eta*H^3*K*v_c/(3*L_K^3) with L_K = 180 km
(the H-sweep fit, plots/tmp/test_dqds_Hexponent.py).  Lambda uses the same
scaling.  The output name carries ".n3" for n = 3.

A seventh argument "pooled" draws every plate configuration on one row of
panels instead of one row each (same markers and colours; ".pooled" in the
output name).

Usage:  python3 plot_forces.varH.py 300.e3 10.e3 10.e3 1.e3 [n [L_km [pooled]]]
"""
import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from functions import nominal_plate_thickness

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.labelsize'] = 7
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6

analysis_depth    = float(sys.argv[1])
analysis_depth_dz = float(sys.argv[2])
ds                = float(sys.argv[3])
dz                = float(sys.argv[4])
n_exp = int(sys.argv[5]) if len(sys.argv) > 5 else 1
if n_exp == 1:
    L_eff = float(sys.argv[6]) if len(sys.argv) > 6 else 1497.0   # km
    scal_label = r'$\eta H K v_c / L_{eff}$   [MPa]'
else:
    L_eff = float(sys.argv[6]) if len(sys.argv) > 6 else 180.0    # km, L_K
    scal_label = r'$\eta H^%d K v_c / (%d L_K^%d)$   [MPa]' % (n_exp, n_exp, n_exp)

tactual_min = 11
tmin = tactual_min - 8

suffix = '.z%s.shear-dz%s.ds%s.prof-dz%skm' % (analysis_depth/1.e3, analysis_depth_dz/1.e3, ds/1.e3, dz/1.e3)
pooled_fig = len(sys.argv) > 7 and sys.argv[7] == 'pooled'
plot_name = 'plots/DP-comparisons/compilations/forces-vs-scaling.varH' + suffix + '.tmin%d' % tmin \
            + ('' if n_exp == 1 else '.n%d' % n_exp) + ('.pooled' if pooled_fig else '')

K_ind = 11; ss_ind = 17; sn_ind = 6; anal_ind = 4; DP_ind = 3; vc_ind = 19; cov_ind = 38

mant_visc  = 2.5e20
cmyr_to_ms = 0.01/(365.25*24*3600)
drho, g    = 50., 9.81

H100 = {'free':    {250: '2D_compositional_subd_lower-res_new_250plates',
                    500: '2D_compositional_subd_lower-res_new2',
                    1000: '2D_compositional_subd_lower-res_new_1000plates'},
        'FixedSP': {250: '2D_compositional_subd_lower-res_new_FixedSP_250plates',
                    500: '2D_compositional_subd_FixedSP_lower-res_new2',
                    1000: '2D_compositional_subd_lower-res_new_FixedSP_1000plates2'},
        'FixedOP': {250: '2D_compositional_subd_lower-res_new_FixedOP_250plates',
                    500: '2D_compositional_subd_FixedOP_lower-res_new',
                    1000: '2D_compositional_subd_lower-res_new_FixedOP_1000plates'}}

def name(bc, eta, H):
    if bc == 'drho71':
        return '2D_compositional_subd_lower-res_new_500plates_H70km_drho71'
    if H == 100:
        return H100[bc][eta]
    return '2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % ('' if bc == 'free' else bc + '_', eta, H)

BCS  = ['free', 'FixedSP', 'FixedOP']
ETAS = [250, 500, 1000]
HS   = [70, 100, 130]
bc_label = {'free': 'free plates', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
colors  = {250: 'peru', 500: 'maroon', 1000: 'black'}
markers = {70: 'v', 100: 'o', 130: 's'}
# Overturned slabs, plotted hollow and left out of the statistics as in Fig. 6.
# The H = 100 entries are the paper's classification (eta' = 1000 free, and
# eta' = 500, 1000 fixed OP).  The two free H = 130 runs are judged from their
# final shapes (rolled over in the lower mantle, like 1000/100) and dip
# histories (shallowing to 51-56 deg with K up to 0.0024-0.0027 /km); their
# eta'*H^3 is 1.1x and 2.2x that of 1000/100.  Fixed-plate H = 70/130 runs are
# classified as they come in (2026-10-08: fixed OP 250/500/1000 at H = 70
# checked, all normal).
# 2026-10-10, from plots/tmp/varH_final_shapes.png: every fixed-SP slab is
# straight; all three fixed-OP H = 130 slabs and fixed-OP 1000/70 fold back
# toward the SP in the lower mantle (1000/70 is the borderline call).
OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130),
              ('FixedOP', 500, 100), ('FixedOP', 1000, 100),
              ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70)}
# the B-matched run: drho = 71.43 so that B = drho*g*H = 49 MPa at H = 70 km
DRHO = {'drho71': 71.43}

def load(bc, eta, H):
    m = name(bc, eta, H)
    f = 'text_files/withT_varH/' + m + suffix + '.txt'
    if not os.path.exists(f):
        print('skipping %s (no withT_varH file)' % m)
        return None
    d = np.loadtxt(f)[tmin:]
    d = d[d[:, cov_ind] == 1]
    Hn = nominal_plate_thickness(m)                       # km
    B = DRHO.get(bc, drho) * g * Hn * 1.e3 * 1e-6         # MPa
    # n = 1: eta*H*K*v_c/L_eff;  n > 1: eta*H^n*K*v_c/(n*L^n)   [MPa]
    scal = (Hn/L_eff)**n_exp / n_exp**(n_exp > 1) * d[:, K_ind] * eta * mant_visc * d[:, vc_ind] * cmyr_to_ms * 1e-6
    dqds = -d[:, sn_ind] / 1.e6                           # MPa, paper sign
    dpmis = -(d[:, DP_ind] - d[:, anal_ind]) / 1.e6       # B_slab - DP, MPa
    DP = d[:, DP_ind] / 1.e6                              # MPa
    Bslab = d[:, anal_ind] / 1.e6                         # drho*g*H*cos(dip), MPa
    full = (d[:, DP_ind] - d[:, anal_ind] + d[:, ss_ind] - d[:, sn_ind]) / 1.e6
    return dict(scal=scal, dqds=dqds, dpmis=dpmis, full=full, B=B, Hn=Hn, DP=DP, Bslab=Bslab)

data = {(bc, e, H): load(bc, e, H) for bc in BCS for e in ETAS for H in HS}
data[('drho71', 500, 70)] = load('drho71', 500, 70)
data = {k: v for k, v in data.items() if v is not None}
bcs_present = [bc for bc in BCS if any(k[0] == bc for k in data)]

# --- statistics
def implied_L(d):
    ok = (d['scal'] > 0) & (d['dqds'] > 0)
    return d['scal'][ok] / d['dqds'][ok] * L_eff, ok.sum()

print('%-8s %-6s %-5s %4s %7s %9s %9s %9s' % ('plates', 'eta', 'H', 'N', 'r', 'L_impl', 'med|full|', 'B'))
for (bc, e, H), d in data.items():
        r = pearsonr(d['scal'], d['dqds'])[0] if len(d['scal']) > 2 else np.nan
        L, n = implied_L(d)
        tag = ' (overturned)' if (bc, e, H) in OVERTURNED else ''
        print('%-8s %-6d %-5d %4d %7.3f %9.0f %9.2f %9.1f%s'
              % (bc, e, H, len(d['scal']), r, np.median(L) if n else np.nan,
                 np.median(np.abs(d['full'])), d['B'], tag))
print()

def pooled(label, keys):
    ds_ = [data[k] for k in keys if k in data and k not in OVERTURNED]
    if not ds_:
        return
    scal = np.concatenate([d['scal'] for d in ds_]); dqds = np.concatenate([d['dqds'] for d in ds_])
    Ls = np.concatenate([implied_L(d)[0] for d in ds_])
    # slope through the origin of dQ/ds on the scaling: 1 if the length scale
    # holds here, and L_eff/slope (n = 1) or L_K/slope^(1/n) (n > 1) is the
    # best-fitting length
    slope = np.sum(scal*dqds) / np.sum(scal**2)
    B = np.median([d['B'] for d in ds_])
    print('  %-18s N=%3d (%d models)  r=%.3f  slope=%.2f (L_fit=%4.0f km)  rms(dQ/ds-scal)/B=%.3f  '
          'implied L_eff median=%5.0f km (N=%d, IQR %.0f-%.0f)'
          % (label, len(scal), len(ds_), pearsonr(scal, dqds)[0], slope, L_eff/slope**(1./n_exp),
             np.sqrt(np.mean((dqds - scal)**2))/B, np.median(Ls), len(Ls),
             np.percentile(Ls, 25), np.percentile(Ls, 75)))

print('per plate configuration and H, non-overturned models pooled:')
for bc in bcs_present:
    for H in HS:
        pooled('%s H=%d' % (bc, H), [(bc, e, H) for e in ETAS])
print('per H, all plate configurations pooled:')
for H in HS:
    pooled('H=%d' % H, [(bc, e, H) for bc in BCS for e in ETAS])
ds_ = [data[k] for k in data if k not in OVERTURNED]
scal = np.concatenate([d['scal'] for d in ds_]); dqds = np.concatenate([d['dqds'] for d in ds_])
lam = np.concatenate([d['scal']/d['B'] for d in ds_]); q = np.concatenate([d['dqds']/d['B'] for d in ds_])
print('  all H  N=%3d  r=%.3f (MPa)  r=%.3f (over B)  slope=%.2f'
      % (len(scal), pearsonr(scal, dqds)[0], pearsonr(lam, q)[0], np.sum(scal*dqds)/np.sum(scal**2)))

# --- figure: one row of three panels per plate configuration, or one row for all
rows = ['all'] if pooled_fig else bcs_present
nrow = len(rows)
fig, axs = plt.subplots(nrow, 3, figsize=(7.0, 2.6*nrow), squeeze=False)

def scatter(ax, x, y, bc, e, H):
    if bc == 'drho71':       # B-matched run: free plates, drawn as diamonds
        ax.scatter(x, y, s=12, color=colors[e], edgecolor='white', linewidth=0.4,
                   marker='D', zorder=4)
    elif (bc, e, H) in OVERTURNED:
        ax.scatter(x, y, s=10, facecolors='none', edgecolor=colors[e], linewidth=0.4,
                   marker=markers[H], alpha=0.7, zorder=3)
    else:
        ax.scatter(x, y, s=10, color=colors[e], edgecolor='black', linewidth=0.25,
                   marker=markers[H], zorder=3)

for (bc, e, H), d in data.items():
    row = 0 if pooled_fig else bcs_present.index('free' if bc == 'drho71' else bc)
    scatter(axs[row, 0], d['Bslab'], d['DP'], bc, e, H)
    scatter(axs[row, 1], d['scal'], d['dqds'], bc, e, H)
    scatter(axs[row, 2], d['scal'], d['dpmis'], bc, e, H)

lims = [(-5, 45, -5, 45), (-5, 40, -10, 30), (-5, 40, -15, 30)]
xlabels = [r'$B_{slab}$   [MPa]', scal_label, scal_label]
ylabels = [r'$\Delta P$   [MPa]', r'$dQ_s/ds$   [MPa]', r'$B_{slab} - \Delta P$   [MPa]']
for row, bc in enumerate(rows):
    for ax, (x0, x1, y0, y1), xl, yl in zip(axs[row], lims, xlabels, ylabels):
        xx = np.linspace(x0, x1, 10)
        ax.plot(xx, xx, color='gray', linestyle='--', linewidth=1, zorder=-1)
        ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
        ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
        ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.set_aspect((x1 - x0)/(y1 - y0), adjustable='box')
    axs[row, 0].set_title(', '.join(bc_label[b] for b in bcs_present) if bc == 'all' else bc_label[bc],
                          loc='left', fontsize=7)

from matplotlib.lines import Line2D
handles = [Line2D([], [], marker=markers[H], color='gray', linestyle='', label='H = %d km' % H) for H in HS] + \
          [Line2D([], [], marker='o', color=colors[e], linestyle='', label=r"$\eta'$ = %d" % e) for e in ETAS]
if ('drho71', 500, 70) in data:
    handles.append(Line2D([], [], marker='D', color=colors[500], linestyle='',
                          label=r'H = 70, $\Delta\rho$ = 71 (B = 49 MPa)'))
axs[0, 0].legend(handles=handles, loc='upper left', fontsize=5, frameon=False)

plt.tight_layout()
plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
print('wrote', plot_name + '.png/.pdf')
