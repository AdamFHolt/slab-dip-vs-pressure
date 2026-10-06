#!/bin/python
"""
Variable plate thickness check (reviewer 2, comment 2): does the dQ/ds scaling
eta*H*K*v_c/L_eff, calibrated on H = 100 km, hold at H = 70 and 130 km?

Nine free-plate models: eta' = 250, 500, 1000 crossed with H = 70, 100, 130 km.
Reads text_files/withT_varH/ (built by make_withT.varH.py), whose H = 100 files
come from extracted_coverage/ and the others from extracted_varH/, so all nine
are from the same version of extract_properties.py.  Timesteps whose midplane
does not resolve the analysis depth (coverage flag, column 38) are dropped, as
in withT_covered/.

  (a) dQ/ds against the scaling eta*H*K*v_c/L_eff  [MPa]
  (b) the same, both divided by B = drho*g*H  (Lambda on the x-axis)
  (c) (B_slab - DP)/B against Lambda

Writes plots/DP-comparisons/compilations/forces-vs-scaling.varH.z<d>...png/pdf
and prints, per model and per H, the correlation and the implied L_eff
(median of eta*H*K*v_c / (dQ/ds) over timesteps where both are positive, the
same sign filter as the L_eff = 1497 km calibration).

Usage:  python3 plot_forces.varH.py 300.e3 10.e3 10.e3 1.e3
"""
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
L_eff = 1497.0  # km

tactual_min = 11
tmin = tactual_min - 8

suffix = '.z%s.shear-dz%s.ds%s.prof-dz%skm' % (analysis_depth/1.e3, analysis_depth_dz/1.e3, ds/1.e3, dz/1.e3)
plot_name = 'plots/DP-comparisons/compilations/forces-vs-scaling.varH' + suffix + '.tmin%d' % tmin

K_ind = 11; ss_ind = 17; sn_ind = 6; anal_ind = 4; DP_ind = 3; vc_ind = 19; cov_ind = 38

mant_visc  = 2.5e20
cmyr_to_ms = 0.01/(365.25*24*3600)
drho, g    = 50., 9.81

def name(eta, H):
    if H == 100:
        return {250: '2D_compositional_subd_lower-res_new_250plates',
                500: '2D_compositional_subd_lower-res_new2',
                1000: '2D_compositional_subd_lower-res_new_1000plates'}[eta]
    return '2D_compositional_subd_lower-res_new_%dplates_H%dkm' % (eta, H)

ETAS = [250, 500, 1000]
HS   = [70, 100, 130]
colors  = {250: 'peru', 500: 'maroon', 1000: 'black'}
markers = {70: 'v', 100: 'o', 130: 's'}
# Overturned slabs, plotted hollow and left out of the statistics as in Fig. 6.
# (1000, 100) is the paper's classification.  The two H = 130 runs are judged
# from their final shapes (rolled over in the lower mantle, like 1000/100) and
# dip histories (shallowing to 51-56 deg with K up to 0.0024-0.0027 /km);
# their eta'*H^3 is 1.1x and 2.2x that of 1000/100.
OVERTURNED = {(1000, 100), (500, 130), (1000, 130)}

def load(eta, H):
    m = name(eta, H)
    d = np.loadtxt('text_files/withT_varH/' + m + suffix + '.txt')[tmin:]
    d = d[d[:, cov_ind] == 1]
    Hn = nominal_plate_thickness(m)                       # km
    B = drho * g * Hn * 1.e3 * 1e-6                       # MPa
    scal = Hn/L_eff * d[:, K_ind] * eta * mant_visc * d[:, vc_ind] * cmyr_to_ms * 1e-6   # MPa
    dqds = -d[:, sn_ind] / 1.e6                           # MPa, paper sign
    dpmis = -(d[:, DP_ind] - d[:, anal_ind]) / 1.e6       # B_slab - DP, MPa
    full = (d[:, DP_ind] - d[:, anal_ind] + d[:, ss_ind] - d[:, sn_ind]) / 1.e6
    return dict(scal=scal, dqds=dqds, dpmis=dpmis, full=full, B=B, Hn=Hn)

data = {(e, H): load(e, H) for e in ETAS for H in HS}

# --- statistics
def implied_L(d):
    ok = (d['scal'] > 0) & (d['dqds'] > 0)
    return d['scal'][ok] / d['dqds'][ok] * L_eff, ok.sum()

print('%-6s %-5s %4s %7s %9s %9s %9s' % ('eta', 'H', 'N', 'r', 'L_impl', 'med|full|', 'B'))
for e in ETAS:
    for H in HS:
        d = data[(e, H)]
        r = pearsonr(d['scal'], d['dqds'])[0] if len(d['scal']) > 2 else np.nan
        L, n = implied_L(d)
        tag = ' (overturned)' if (e, H) in OVERTURNED else ''
        print('%-6d %-5d %4d %7.3f %9.0f %9.2f %9.1f%s'
              % (e, H, len(d['scal']), r, np.median(L) if n else np.nan,
                 np.median(np.abs(d['full'])), d['B'], tag))
print()
print('per H, non-overturned models pooled:')
for H in HS:
    ds_ = [data[(e, H)] for e in ETAS if (e, H) not in OVERTURNED]
    scal = np.concatenate([d['scal'] for d in ds_]); dqds = np.concatenate([d['dqds'] for d in ds_])
    Ls = np.concatenate([implied_L(d)[0] for d in ds_])
    # slope through the origin of dQ/ds on the scaling: 1 if L_eff = 1497 km
    # holds at this H, and L_eff/slope is the best-fitting L_eff
    slope = np.sum(scal*dqds) / np.sum(scal**2)
    B = ds_[0]['B']
    print('  H=%3d  N=%3d  r=%.3f  slope=%.2f (L_fit=%4.0f km)  rms(dQ/ds-scal)/B=%.3f  '
          'implied L_eff median=%5.0f km (N=%d, IQR %.0f-%.0f)'
          % (H, len(scal), pearsonr(scal, dqds)[0], slope, L_eff/slope,
             np.sqrt(np.mean((dqds - scal)**2))/B, np.median(Ls), len(Ls),
             np.percentile(Ls, 25), np.percentile(Ls, 75)))
ds_ = [data[k] for k in data if k not in OVERTURNED]
scal = np.concatenate([d['scal'] for d in ds_]); dqds = np.concatenate([d['dqds'] for d in ds_])
lam = np.concatenate([d['scal']/d['B'] for d in ds_]); q = np.concatenate([d['dqds']/d['B'] for d in ds_])
print('  all H  N=%3d  r=%.3f (MPa)  r=%.3f (over B)  slope=%.2f'
      % (len(scal), pearsonr(scal, dqds)[0], pearsonr(lam, q)[0], np.sum(scal*dqds)/np.sum(scal**2)))

# --- figure
fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6))

def scatter(ax, x, y, e, H):
    if (e, H) in OVERTURNED:
        ax.scatter(x, y, s=10, facecolors='none', edgecolor=colors[e], linewidth=0.4,
                   marker=markers[H], alpha=0.7, zorder=3)
    else:
        ax.scatter(x, y, s=10, color=colors[e], edgecolor='black', linewidth=0.25,
                   marker=markers[H], zorder=3)

for (e, H), d in data.items():
    scatter(axs[0], d['scal'], d['dqds'], e, H)
    scatter(axs[1], d['scal']/d['B'], d['dqds']/d['B'], e, H)
    scatter(axs[2], d['scal']/d['B'], d['dpmis']/d['B'], e, H)

lims = [(-5, 40, -5, 40), (-0.1, 0.8, -0.1, 0.8), (-0.1, 0.8, -0.25, 0.4)]
xlabels = [r'$\eta H K v_c / L_{eff}$   [MPa]', r'$\Lambda$', r'$\Lambda$']
ylabels = [r'$dQ_s/ds$   [MPa]', r'$(dQ_s/ds)\,/\,B$', r'$(B_{slab} - \Delta P)\,/\,B$']
for ax, (x0, x1, y0, y1), xl, yl in zip(axs, lims, xlabels, ylabels):
    xx = np.linspace(x0, x1, 10)
    ax.plot(xx, xx, color='gray', linestyle='--', linewidth=1, zorder=-1)
    ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
    ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_aspect((x1 - x0)/(y1 - y0), adjustable='box')

from matplotlib.lines import Line2D
handles = [Line2D([], [], marker=markers[H], color='gray', linestyle='', label='H = %d km' % H) for H in HS] + \
          [Line2D([], [], marker='o', color=colors[e], linestyle='', label=r"$\eta'$ = %d" % e) for e in ETAS]
axs[0].legend(handles=handles, loc='upper left', fontsize=5, frameon=False)

plt.tight_layout()
plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
print('wrote', plot_name + '.png/.pdf')
