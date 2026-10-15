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

Usage:  python3 plot_forces.varH.py 300.e3 10.e3 10.e3 1.e3 [n [L_km [pooled] [nodrho71] [alleta] [closure|closureK]]]
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
pooled_fig = 'pooled' in sys.argv[7:]
# the B-matched drho71 run is a supplement control (it shows that thickness, not
# buoyancy, sets dQ/ds); "nodrho71" leaves it out for the main-text figure
with_drho71 = 'nodrho71' not in sys.argv[7:]
# "alleta" adds the paper's eta' = 50 and 375 models, which exist only at H = 100
all_eta = 'alleta' in sys.argv[7:]
# "closure": Fig. 3 instead -- DP + sigma_slab vs B_slab, one panel per plate
# thickness.  "closureK": its supplement companion, the residual
# DP + sigma_slab - B_slab against |K| in the same three panels.  Every output.
closure_fig  = 'closure' in sys.argv[7:]
closureK_fig = 'closureK' in sys.argv[7:]
# "fig5": Fig. 5 -- DP vs B_slab, dQ/ds vs the scaling and B_slab - DP vs the
# scaling (n and L from arguments 5-6) side by side, every model of every H in
# the same panels, in the style of the submitted figure (plot_forces.no-ot.py)
fig5_fig     = 'fig5' in sys.argv[7:]
# the closure figures (Fig. 3 and its supplement versions) go to manuscript_figs/
out_dir = 'plots/DP-comparisons/compilations/' + ('manuscript_figs/' if (closure_fig or closureK_fig or fig5_fig) else '')
os.makedirs(out_dir, exist_ok=True)
plot_name = out_dir + 'forces-vs-scaling.varH' + suffix + '.tmin%d' % tmin \
            + ('' if n_exp == 1 else '.n%d' % n_exp) + ('.pooled' if pooled_fig else '') + ('' if with_drho71 else '.nodrho71') + ('.alleta' if all_eta else '') + ('.closure' if closure_fig else '') + ('.closureK' if closureK_fig else '') + ('.fig5' if fig5_fig else '')

K_ind = 11; ss_ind = 17; sn_ind = 6; anal_ind = 4; DP_ind = 3; vc_ind = 19; cov_ind = 38

mant_visc  = 2.5e20
cmyr_to_ms = 0.01/(365.25*24*3600)
drho, g    = 50., 9.81

H100 = {'free':    {50: '2D_compositional_subd_lower-res_new_50plates',
                    250: '2D_compositional_subd_lower-res_new_250plates',
                    375: '2D_compositional_subd_lower-res_new_375plates',
                    500: '2D_compositional_subd_lower-res_new2',
                    1000: '2D_compositional_subd_lower-res_new_1000plates'},
        'FixedSP': {50: '2D_compositional_subd_lower-res_new_FixedSP_50plates',
                    250: '2D_compositional_subd_lower-res_new_FixedSP_250plates',
                    375: '2D_compositional_subd_lower-res_new_FixedSP_375plates',
                    500: '2D_compositional_subd_FixedSP_lower-res_new2',
                    1000: '2D_compositional_subd_lower-res_new_FixedSP_1000plates2'},
        'FixedOP': {50: '2D_compositional_subd_lower-res_new_FixedOP_50plates',
                    250: '2D_compositional_subd_lower-res_new_FixedOP_250plates',
                    375: '2D_compositional_subd_lower-res_new_FixedOP_375plates',
                    500: '2D_compositional_subd_FixedOP_lower-res_new',
                    1000: '2D_compositional_subd_lower-res_new_FixedOP_1000plates'}}

def name(bc, eta, H):
    if bc == 'drho71':
        return '2D_compositional_subd_lower-res_new_500plates_H70km_drho71'
    if H == 100:
        return H100[bc][eta]
    return '2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % ('' if bc == 'free' else bc + '_', eta, H)

BCS  = ['free', 'FixedSP', 'FixedOP']
ETAS = [50, 250, 375, 500, 1000] if all_eta else [250, 500, 1000]
HS   = [70, 100, 130]
bc_label = {'free': 'free plates', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
colors  = {50: 'tan', 250: 'peru', 375: 'firebrick', 500: 'maroon', 1000: 'black'}
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
              ('FixedOP', 1000, 70),
              ('FixedOP', 375, 100)}       # as in the paper
# the B-matched run: drho = 71.43 so that B = drho*g*H = 49 MPa at H = 70 km
DRHO = {'drho71': 71.43}
# last output used, where it is earlier than the extraction window's end: the
# eta' = 50 fixed-SP slab necks at 250 km on its final output (thickness 66 ->
# 44 km, dQ/ds jumping to 6.5 MPa), so that output is left out (2026-10-12)
LAST_OUTPUT = {('FixedSP', 50, 100): 22}

def load(bc, eta, H):
    m = name(bc, eta, H)
    f = 'text_files/withT_varH/' + m + suffix + '.txt'
    if not os.path.exists(f):
        if not (eta in (50, 375) and H != 100):     # those were only ever run at H = 100
            print('skipping %s (no withT_varH file)' % m)
        return None
    d = np.loadtxt(f)[tmin:]
    d = d[d[:, cov_ind] == 1]
    d = d[d[:, 0] <= LAST_OUTPUT.get((bc, eta, H), 1e9)]
    Hn = nominal_plate_thickness(m)                       # km
    B = DRHO.get(bc, drho) * g * Hn * 1.e3 * 1e-6         # MPa
    # n = 1: eta*H*K*v_c/L_eff;  n > 1: eta*H^n*K*v_c/(n*L^n)   [MPa]
    scal = (Hn/L_eff)**n_exp / n_exp**(n_exp > 1) * d[:, K_ind] * eta * mant_visc * d[:, vc_ind] * cmyr_to_ms * 1e-6
    dqds = -d[:, sn_ind] / 1.e6                           # MPa, paper sign
    dpmis = -(d[:, DP_ind] - d[:, anal_ind]) / 1.e6       # B_slab - DP, MPa
    DP = d[:, DP_ind] / 1.e6                              # MPa
    Bslab = d[:, anal_ind] / 1.e6                         # drho*g*H*cos(dip), MPa
    full = (d[:, DP_ind] - d[:, anal_ind] + d[:, ss_ind] - d[:, sn_ind]) / 1.e6
    K = np.abs(d[:, K_ind]) * 1.e6                        # |K|, 1e-3/km, as in Fig. 3
    return dict(scal=scal, dqds=dqds, dpmis=dpmis, full=full, B=B, Hn=Hn, DP=DP, Bslab=Bslab, K=K)

data = {(bc, e, H): load(bc, e, H) for bc in BCS for e in ETAS for H in HS}
if with_drho71:
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

# --- closure figures (Fig. 3 and its supplement companion): one panel per H.
# Colour is eta', marker the plate configuration, hollow = overturned.  The median
# |residual| of the non-overturned models and any points outside the axes are
# printed, not drawn (they belong in the caption).  Text is larger than the
# three-row figures so the single row reads at the same size in print.
if closure_fig or closureK_fig:
    mpl.rcParams.update({'axes.labelsize': 9, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
    from matplotlib.lines import Line2D
    bc_marker = {'FixedSP': 'v', 'free': 'o', 'FixedOP': '^', 'drho71': 'D'}
    if closure_fig:
        # Fig. 3 (300 km) starts at -2.5 MPa; the supplement depths keep -5, which
        # their lowest points need
        lo = -2.5 if analysis_depth == 300.e3 else -5
        (x0, x1, y0, y1) = (lo, 40, lo, 40)
        xlab, ylab = r'$B_{slab}$   [MPa]', r'$\Delta P + \sigma_{slab}$   [MPa]'
    else:
        (x0, x1, y0, y1) = (-0.2, 3.2, -10, 10)
        xlab, ylab = r'$|K|$   [10$^{-3}$ km$^{-1}$]', r'$\Delta P + \sigma_{slab} - B_{slab}$   [MPa]'
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6), sharey=True)
    for ax, H, lab in zip(axs, HS, 'abc'):
        keys = [k for k in data if k[2] == H]
        n_off = {True: 0, False: 0}     # keyed by 'overturned'
        for (bc, e, Hk) in keys:
            d = data[(bc, e, Hk)]
            x = d['Bslab'] if closure_fig else d['K']
            y = d['full'] + d['Bslab'] if closure_fig else d['full']
            n_off[(bc, e, Hk) in OVERTURNED] += np.sum((x < x0) | (x > x1) | (y < y0) | (y > y1))
            m = bc_marker[bc]
            if (bc, e, Hk) in OVERTURNED:
                ax.scatter(x, y, s=10, facecolors='none', edgecolor=colors[e], linewidth=0.4,
                           marker=m, alpha=0.7, zorder=3)
            else:
                ax.scatter(x, y, s=10, color=colors[e], edgecolor='black' if bc != 'drho71' else 'white',
                           linewidth=0.25, marker=m, zorder=3)
        res = np.concatenate([data[k]['full'] for k in keys if k not in OVERTURNED])
        if closure_fig:
            ax.plot([x0, x1], [x0, x1], color='gray', linestyle='--', linewidth=1, zorder=1)
            ax.set_xticks([0, 10, 20, 30, 40]); ax.set_yticks([0, 10, 20, 30, 40])
        else:
            ax.set_yticks([-10, -5, 0, 5, 10])
        ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
        ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
        ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect((x1 - x0)/(y1 - y0), adjustable='box')
        ax.set_title('H = %d km' % H, fontsize=9)
        # panel label on the title's baseline, flush with the axes' left edge
        ax.set_title(lab + ')', loc='left', fontsize=9)
        ax.set_xlabel(xlab)
        print('closure H=%d: N=%d normal, median |res| %.2f MPa; off scale: %d normal, %d overturned'
              % (H, len(res), np.median(np.abs(res)), n_off[False], n_off[True]))
    axs[0].set_ylabel(ylab)
    # two columns: eta' colours on the left, plate-configuration markers and the
    # overturned symbol on the right (a blank pads the left column to four rows)
    blank = Line2D([], [], linestyle='', label=' ')
    eta_h = [Line2D([], [], marker='o', color=colors[e], markeredgecolor='black', markeredgewidth=0.25,
                    linestyle='', markersize=4, label="\u03b7' = %d" % e) for e in ETAS]
    # plain-text eta (Myriad has the glyph, not the prime) so both columns share one row height
    bc_h = [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
                   linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
           [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray',
                   linestyle='', markersize=4, label='overturned')]
    if ('drho71', 500, 70) in data:
        eta_h.append(Line2D([], [], marker='D', color=colors[500], linestyle='', markersize=4,
                            label=r'$\Delta\rho$ = 71'))
    eta_h += [blank] * (len(bc_h) - len(eta_h))
    # Fig. 3: the empty upper left; the residual panels: the empty lower half
    leg = axs[0].legend(handles=eta_h + bc_h, ncol=2, loc='upper left' if closure_fig else 'lower left',
                        fontsize=6.5, frameon=True, fancybox=False, edgecolor='black', facecolor='white',
                        framealpha=1, borderpad=0.5, labelspacing=0.35, handletextpad=0.1,
                        handlelength=1.2, columnspacing=0.8, borderaxespad=0.6)
    leg.get_frame().set_linewidth(0.5)
    plt.tight_layout()
    plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(plot_name + '.pdf', bbox_inches='tight')
    print('wrote', plot_name + '.png/.pdf')
    sys.exit(0)

# --- Fig. 5, all models pooled.  Marker styles as plot_forces.no-ot.py and
# functions_plotting.py (s = 10, thin black edge, overturned hollow); colour
# eta', marker the plate configuration; legend as in Fig. 3.  The Lambda top
# axes of the submitted figure are left off: Lambda = scaling/B and B differs
# between the three H.
if fig5_fig:
    from matplotlib.lines import Line2D
    mpl.rcParams.update({'axes.labelsize': 9, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
    bc_marker = {'FixedSP': 'v', 'free': 'o', 'FixedOP': '^'}
    fig, (axP, axQ, axM) = plt.subplots(1, 3, figsize=(7.0, 2.6))
    for (bc, e, H), d in data.items():
        m = bc_marker[bc]
        for ax, x, y in ((axP, d['Bslab'], d['DP']), (axQ, d['scal'], d['dqds']), (axM, d['scal'], d['dpmis'])):
            if (bc, e, H) in OVERTURNED:
                ax.scatter(x, y, s=10, facecolors='none', edgecolor=colors[e], marker=m,
                           linewidth=0.4, alpha=0.7, zorder=2)
            else:
                ax.scatter(x, y, s=10, color=colors[e], edgecolor='black', marker=m,
                           linewidth=0.25, zorder=4 if (ax is axP and bc == 'FixedSP') else 3)
    for ax, (x0, x1, y0, y1), xl, yl, lab in (
            (axP, (-3, 40, -3, 40), r'$B_{slab}$   [MPa]', r'$\Delta P$   [MPa]', 'a'),
            (axQ, (-5, 40, -10, 30), scal_label, r'$-dQ_{s}/ds$   [MPa]', 'b'),
            (axM, (-5, 40, -10, 30), scal_label, r'$-(\Delta P - B_{slab})$   [MPa]', 'c')):
        if ax is axP:
            ax.plot([x0, x1], [x0, x1], color='black', linewidth=1, zorder=2)
            ax.set_yticks([0, 10, 20, 30, 40])
        else:
            ax.plot([x0, x1], [x0, x1], color='gray', linestyle='--', linewidth=1, zorder=1)
        ax.set_xticks([0, 10, 20, 30, 40])
        ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
        ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
        ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect((x1 - x0)/(y1 - y0), adjustable='box')
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.text(0.04, 0.96, lab + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10)
    # legend: one boxed row above all three panels (inside any panel it would
    # cover points); eta' colours, then plate configuration and overturned
    eta_h = [Line2D([], [], marker='o', color=colors[e], markeredgecolor='black', markeredgewidth=0.25,
                    linestyle='', markersize=4, label="η' = %d" % e) for e in ETAS]
    bc_h = [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
                   linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
           [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray',
                   linestyle='', markersize=4, label='overturned')]
    plt.tight_layout()
    leg = fig.legend(handles=eta_h + bc_h, ncol=len(eta_h + bc_h), loc='lower center',
                     bbox_to_anchor=(0.5, 0.91), fontsize=6.5, frameon=True, fancybox=False,
                     edgecolor='black', facecolor='white', framealpha=1, borderpad=0.5,
                     handletextpad=0.1, handlelength=1.2, columnspacing=1.0)
    leg.get_frame().set_linewidth(0.5)
    plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(plot_name + '.pdf', bbox_inches='tight')
    print('wrote', plot_name + '.png/.pdf')
    sys.exit(0)

# --- figure: one row of three panels per plate configuration, or one row for all
rows = ['all'] if pooled_fig else bcs_present
nrow = len(rows)
ncol = 3
fig, axs = plt.subplots(nrow, ncol, figsize=(7.0/3*ncol, 2.6*nrow), squeeze=False)

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
one_to_one = [True, True, True]
for row, bc in enumerate(rows):
    for ax, (x0, x1, y0, y1), xl, yl, oto in zip(axs[row], lims, xlabels, ylabels, one_to_one):
        xx = np.linspace(x0, x1, 10)
        if oto:
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
