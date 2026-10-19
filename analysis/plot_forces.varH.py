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
# "fig5lam": the same with panels b) and c) in Lambda space -- x = scaling/B,
# y = (-dQ/ds)/B and (B_slab - DP)/B, each model over its own nominal
# B = drho*g*H (the Earth maps' definition), so every H shares one set of axes
fig5lam_fig  = 'fig5lam' in sys.argv[7:]
fig5_fig     = fig5_fig or fig5lam_fig
# "means" (with fig5/fig5lam): panels b) and c) show every output pale in the
# background and each model's time mean on top, with +-1 sigma whiskers; panel
# a) keeps every output (DP tracks B_slab in time, which is the point there)
means_fig    = fig5_fig and 'means' in sys.argv[7:]
# "inset" (with means): a small inset in panels b) and c) with every model
# mean, overturned (hollow) included, on axes wide enough for them
inset_fig    = means_fig and 'inset' in sys.argv[7:]
# "lamd" (with fig5, not fig5lam): panel c) becomes the residual fraction
# (B_slab - DP)/B against Lambda = scaling/B (nominal B = drho*g*H), with the
# model-mean inset like b)'s -- the quantity the Earth maps use.  The
# absolute residual panel is dropped (it is this one times each model's B)
lamd_fig     = fig5_fig and not fig5lam_fig and 'lamd' in sys.argv[7:]
# "Hsize" (with means): the model means' marker size encodes the plate
# thickness, H = 70 / 100 / 130 km small / medium / large (insets too)
Hsize_fig    = means_fig and 'Hsize' in sys.argv[7:]
MEAN_S       = {70: 16, 100: 30, 130: 50} if Hsize_fig else {70: 30, 100: 30, 130: 30}
INSET_S      = {70: 5, 100: 9, 130: 15} if Hsize_fig else {70: 9, 100: 9, 130: 9}
# "strip" (with lamd, one row): one boxed row like the default, each group
# introduced by a bold header ("viscosity ratio η':" 50 250 ...)
strip_fig    = lamd_fig and 'strip' in sys.argv[7:]
# "lamsupp" (with fig5): supplement figure instead -- Fig. 5a alone, every output
# coloured by |Lambda| (nominal B) rather than eta', so the outputs that leave the
# 1:1 line are seen to be the high-Lambda ones (Adam 2026-10-19: a one-panel
# supplement, not a fourth Fig. 5 panel)
lamsupp_fig  = fig5_fig and 'lamsupp' in sys.argv[7:]
# the manuscript figures go to manuscript_figs/: Fig. 3 and its supplement
# versions (closure, closureK) and Fig. 5 (fig5 means inset lamd Hsize strip);
# the other Fig. 5 variants to compilations/
manuscript_fig = closure_fig or closureK_fig or (fig5_fig and not fig5lam_fig and inset_fig and lamd_fig and Hsize_fig and strip_fig) or lamsupp_fig
out_dir = 'plots/DP-comparisons/compilations/' + ('manuscript_figs/' if manuscript_fig else '')
os.makedirs(out_dir, exist_ok=True)
plot_name = out_dir + 'forces-vs-scaling.varH' + suffix + '.tmin%d' % tmin \
            + ('' if n_exp == 1 else '.n%d' % n_exp) + ('.L%d' % L_eff if fig5_fig else '') + ('.pooled' if pooled_fig else '') + ('' if with_drho71 else '.nodrho71') + ('.alleta' if all_eta else '') + ('.closure' if closure_fig else '') + ('.closureK' if closureK_fig else '') + (('.fig5lam' if fig5lam_fig else '.fig5') if fig5_fig else '') + ('.means' if means_fig else '') + ('.inset' if inset_fig else '') + ('.lamd' if lamd_fig else '') + ('.Hsize' if Hsize_fig else '') + ('.strip' if strip_fig else '')

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
    LAM_LIMS, LAM_TICKS = (-0.1, 0.8, -0.2, 0.6), [0, 0.2, 0.4, 0.6, 0.8]
    if fig5lam_fig or lamd_fig:
        # residual fraction (B_slab - DP)/B in bins of Lambda, normal runs only:
        # where it stays under ~0.1 is the threshold the Earth maps can use
        ok = [d for k, d in data.items() if k not in OVERTURNED]
        lam = np.concatenate([d['scal']/d['B'] for d in ok])
        res = np.concatenate([d['dpmis']/d['B'] for d in ok])
        Hs  = np.concatenate([np.full(len(d['scal']), d['Hn']) for d in ok])
        print('Lambda range %.3f to %.3f; (B_slab - DP)/B against Lambda, normal runs:' % (lam.min(), lam.max()))
        print('  %-13s %4s %7s %7s %7s %7s  %s' % ('Lambda bin', 'N', 'median', 'p75', 'p90', 'f<0.1', 'N by H 70/100/130'))
        for a, b in ((-1, 0), (0, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 2)):
            m = (lam >= a) & (lam < b)
            if m.sum():
                print('  [%5.2f,%5.2f) %4d %7.3f %7.3f %7.3f %7.2f  %s' % (a, b, m.sum(), np.median(res[m]),
                      np.percentile(res[m], 75), np.percentile(res[m], 90), np.mean(res[m] < 0.1),
                      '/'.join(str(np.sum(m & np.isclose(Hs, h, atol=5))) for h in HS)))
        print('  r(Lambda, dQ/ds/B) = %.3f, r(Lambda, residual/B) = %.3f'
              % (pearsonr(lam, np.concatenate([d['dqds']/d['B'] for d in ok]))[0], pearsonr(lam, res)[0]))
    if means_fig:
        # statistics of the model means, normal runs, as drawn (over B in Lambda space)
        ok = [d for k, d in data.items() if k not in OVERTURNED]
        nbs = [d['B'] if fig5lam_fig else 1. for d in ok]
        xs = np.array([np.mean(d['scal'])/b for d, b in zip(ok, nbs)])
        for lab_, key in (('b) dQ/ds', 'dqds'), ('c) B_slab - DP', 'dpmis')):
            ys = np.array([np.mean(d[key])/b for d, b in zip(ok, nbs)])
            print('model means %-15s N=%d  r=%.3f  slope=%.2f  rms=%.2f' % (lab_, len(xs), pearsonr(xs, ys)[0],
                  np.sum(xs*ys)/np.sum(xs**2), np.sqrt(np.mean((ys - xs)**2))))
    if lamsupp_fig:
        # Fig. 5a alone, every output coloured by |Lambda| (the sign is K's; the
        # size is what matters), |Lambda| linear 0-0.5 (above saturates: overturned
        # runs); marker = plate configuration, overturned hollow with a coloured
        # edge.  Colorbar and the plate legend sit in the empty upper left.
        from cmcrameri import cm as cmc
        fig, ax = plt.subplots(figsize=(3.4, 3.2))
        normC = mpl.colors.Normalize(vmin=0, vmax=0.5)
        lam_all = np.abs(np.concatenate([d['scal']/d['B'] for d in data.values()]))
        print('lamsupp: |Lambda| median %.3f, %d of %d outputs > 0.5 (all overturned: %s)'
              % (np.median(lam_all), (lam_all > 0.5).sum(), len(lam_all),
                 all(np.all(np.abs(d['scal']/d['B']) <= 0.5) for k, d in data.items() if k not in OVERTURNED)))
        for (bc, e, H), d in data.items():
            lam = np.abs(d['scal']/d['B'])
            if (bc, e, H) in OVERTURNED:
                ax.scatter(d['Bslab'], d['DP'], s=10, facecolors='none', edgecolor=cmc.navia_r(normC(lam)),
                           marker=bc_marker[bc], linewidth=0.5, alpha=0.8, zorder=2)
            else:
                sc = ax.scatter(d['Bslab'], d['DP'], s=10, c=lam, cmap=cmc.navia_r, norm=normC,
                                edgecolor='black', marker=bc_marker[bc], linewidth=0.25,
                                zorder=4 if bc == 'FixedSP' else 3)
        x0, x1 = -3, 40
        ax.plot([x0, x1], [x0, x1], color='black', linewidth=1, zorder=2)
        ax.set_xticks([0, 10, 20, 30, 40]); ax.set_yticks([0, 10, 20, 30, 40])
        ax.axhline(0, color='lightgray', linewidth=1, zorder=0); ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
        ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
        ax.set_xlim(x0, x1); ax.set_ylim(x0, x1); ax.set_aspect(1, adjustable='box')
        ax.set_xlabel(r'$B_{slab}$   [MPa]'); ax.set_ylabel(r'$\Delta P$   [MPa]')
        # white backing under the colorbar block (hides the gridline behind it)
        ax.add_patch(mpl.patches.Rectangle((0.05, 0.70), 0.50, 0.24, transform=ax.transAxes,
                                           facecolor='white', edgecolor='none', zorder=9))
        cax = ax.inset_axes([0.10, 0.79, 0.36, 0.04]); cax.set_zorder(10)
        cb = fig.colorbar(sc, cax=cax, orientation='horizontal', extend='max')
        cb.set_ticks([0, 0.5]); cb.set_ticklabels(['0', '0.5'])
        cb.ax.xaxis.set_minor_locator(mpl.ticker.FixedLocator([0.1, 0.2, 0.3, 0.4]))
        cb.ax.tick_params(labelsize=6, length=2.5, pad=1); cb.ax.tick_params(which='minor', length=1.5)
        cb.ax.set_title(r'$|\Lambda|$', fontsize=7, pad=4)
        # plate-configuration legend below the colorbar, grey markers as in Fig. 5
        bc_h = [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
                       linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
               [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray', linestyle='',
                       markersize=4, label='overturned')]
        leg = ax.legend(handles=bc_h, loc='upper left', bbox_to_anchor=(0.03, 0.70), fontsize=6.5, frameon=False,
                        handletextpad=0.1, handlelength=1.0, labelspacing=0.3, borderaxespad=0)
        leg.set_zorder(10)
        plot_name = out_dir + 'fig5a-lambda_supp' + suffix + '.tmin%d.n%d.L%d' % (tmin, n_exp, L_eff) \
                    + ('' if with_drho71 else '.nodrho71') + ('.alleta' if all_eta else '')
        plt.savefig(plot_name + '.png', dpi=600, bbox_inches='tight')
        plt.savefig(plot_name + '.pdf', bbox_inches='tight')
        print('wrote', plot_name + '.png/.pdf')
        sys.exit(0)
    if lamd_fig:
        fig, (axP, axQ, axL) = plt.subplots(1, 3, figsize=(7.0, 2.6))
        axM = None
    else:
        fig, (axP, axQ, axM) = plt.subplots(1, 3, figsize=(7.0, 2.6))
        axL = None
    # extents of b), c) (as set below), for the model means drawn unclipped
    lims_m = LAM_LIMS if fig5lam_fig else (-4, 16, -4, 16)
    mean_lims = {axQ: lims_m, axM: lims_m, axL: (-0.1, 0.3, -0.1, 0.3)}
    for (bc, e, H), d in data.items():
        m = bc_marker[bc]
        nb = d['B'] if fig5lam_fig else 1.     # Lambda space: divide by the nominal B
        panels_xy = ((axP, d['Bslab'], d['DP']), (axQ, d['scal']/nb, d['dqds']/nb), (axM, d['scal']/nb, d['dpmis']/nb))
        if lamd_fig:
            panels_xy = panels_xy[:2] + ((axL, d['scal']/d['B'], d['dpmis']/d['B']),)
        for ax, x, y in panels_xy:
            if means_fig and ax is not axP:
                ot = (bc, e, H) in OVERTURNED
                # every output, pale and small, behind the model means; overturned
                # hollow, so the legend's overturned key still applies
                if ot:
                    ax.scatter(x, y, s=5, facecolors='none', edgecolor=colors[e], marker=m, linewidth=0.3,
                               alpha=0.35, zorder=2)
                else:
                    ax.scatter(x, y, s=4, color=colors[e], marker=m, linewidth=0, alpha=0.45, zorder=2)
                # model mean with +-1 sigma whiskers; overturned hollow and behind
                # the normal means.  Markers may run over the axes edge (clip off),
                # so only means inside the axes are drawn; the rest are in the inset
                xm, ym = np.mean(x), np.mean(y)
                x0_, x1_, y0_, y1_ = mean_lims[ax]
                if not (x0_ <= xm <= x1_ and y0_ <= ym <= y1_):
                    continue
                ax.errorbar(xm, ym, xerr=np.std(x), yerr=np.std(y), fmt='none', ecolor=colors[e],
                            elinewidth=0.4, alpha=0.5, capsize=0, zorder=3)
                if ot:
                    ax.scatter(xm, ym, s=MEAN_S[H], facecolors='none', edgecolor=colors[e], marker=m,
                               linewidth=0.75, zorder=4, clip_on=False)
                else:
                    ax.scatter(xm, ym, s=MEAN_S[H], color=colors[e], edgecolor='black', marker=m,
                               linewidth=0.4, zorder=5, clip_on=False)
            elif (bc, e, H) in OVERTURNED:
                ax.scatter(x, y, s=10, facecolors='none', edgecolor=colors[e], marker=m,
                           linewidth=0.4, alpha=0.7, zorder=2)
            else:
                ax.scatter(x, y, s=10, color=colors[e], edgecolor='black', marker=m,
                           linewidth=0.25, zorder=4 if (ax is axP and bc == 'FixedSP') else 3)
    if fig5lam_fig:
        panels_bc = ((axQ, LAM_LIMS, r'$\Lambda$', r'$-(dQ_{s}/ds)\,/\,B$', 'b'),
                     (axM, LAM_LIMS, r'$\Lambda$', r'$-(\Delta P - B_{slab})\,/\,B$', 'c'))
    else:
        # "means" zooms on the normal model means (all within 0-10 MPa)
        lims_bc = (-4, 16, -4, 16) if means_fig else (-5, 40, -10, 30)
        panels_bc = ((axQ, lims_bc, scal_label, r'$-dQ_{s}/ds$   [MPa]', 'b'),
                     (axM, lims_bc, scal_label, r'$-(\Delta P - B_{slab})$   [MPa]', 'c'))
    if lamd_fig:
        panels_bc = panels_bc[:1] + ((axL, (-0.1, 0.3, -0.1, 0.3), r'$\Lambda$', r'$-(\Delta P - B_{slab})\,/\,B$', 'c'),)
    for ax, (x0, x1, y0, y1), xl, yl, lab in (
            ((axP, (-3, 40, -3, 40), r'$B_{slab}$   [MPa]', r'$\Delta P$   [MPa]', 'a'),) + panels_bc):
        if ax is axP:
            ax.plot([x0, x1], [x0, x1], color='black', linewidth=1, zorder=2)
            ax.set_yticks([0, 10, 20, 30, 40])
        else:
            ax.plot([x0, x1], [x0, x1], color='gray', linestyle='--', linewidth=1, zorder=1)
        ax.set_xticks(LAM_TICKS if (fig5lam_fig and ax is not axP) else
                      [-4, 0, 4, 8, 12, 16] if (means_fig and ax is not axP) else [0, 10, 20, 30, 40])
        if means_fig and not fig5lam_fig and ax is not axP:
            ax.set_yticks([-4, 0, 4, 8, 12, 16])
        if ax is axL:
            ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3]); ax.set_yticks([-0.1, 0, 0.1, 0.2, 0.3])
        ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
        ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
        ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        ax.set_aspect((x1 - x0)/(y1 - y0), adjustable='box')
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        # white backing so the label does not sit on the x = 0 line
        ax.text(0.04, 0.96, lab + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10,
                bbox=dict(facecolor='white', edgecolor='none', pad=1))
    if inset_fig:
        # lower right of b) and c), which the normal means leave empty (the inset
        # may run over the right axis): model means only, no background, the main
        # panel's extent outlined
        # (axes, quantity, divide by B, main extent, inset extent and ticks)
        insets = [(axQ, 'dqds', False, lims_bc, (-5, 50, -5, 30), [0, 20, 40], [0, 10, 20, 30])]
        if lamd_fig:
            insets.append((axL, 'dpmis', True, (-0.1, 0.3, -0.1, 0.3), (-0.1, 0.8, -0.1, 0.5),
                           [0, 0.4, 0.8], [0, 0.2, 0.4]))
        else:
            insets.append((axM, 'dpmis', False, lims_bc, (-5, 50, -5, 30), [0, 20, 40], [0, 10, 20, 30]))
        for axO, key, overB, main_lims, inset_lims, ixt, iyt in insets:
            axI = axO.inset_axes([0.62, 0.14, 0.41, 0.36])
            # white backing under the inset and its tick labels, so the pale
            # background points do not run into them
            axO.add_patch(mpl.patches.Rectangle((0.535, 0.045), 0.52, 0.475, transform=axO.transAxes,
                                                facecolor='white', edgecolor='none', clip_on=False, zorder=9))
            for sp in axO.spines.values():       # axes lines over the backing, under the inset
                sp.set_zorder(9.5)
            # short white gaps in the main right axis just above and below the
            # inset, so the axis line does not run straight into its edges
            # (c)'s Lambda inset has its 0.8 tick label right on the main axis
            # line, so its lower gap reaches down past the label)
            for y_gap, h_gap in ((0.50, 0.02), (0.07, 0.07) if overB else (0.12, 0.02)):
                axO.add_patch(mpl.patches.Rectangle((0.95, y_gap), 0.10, h_gap, transform=axO.transAxes,
                                                    facecolor='white', edgecolor='none', clip_on=False, zorder=9.8))
            (ix0, ix1, iy0, iy1) = inset_lims
            for (bc, e, H), d in data.items():
                nb = d['B'] if overB else 1.
                m, xm, ym = bc_marker[bc], np.mean(d['scal'])/nb, np.mean(d[key])/nb
                if (bc, e, H) in OVERTURNED:
                    axI.scatter(xm, ym, s=INSET_S[H], facecolors='white', edgecolor=colors[e], marker=m, linewidth=0.6, zorder=4)
                else:
                    axI.scatter(xm, ym, s=INSET_S[H], color=colors[e], edgecolor='black', marker=m, linewidth=0.25, zorder=5)
            axI.plot([ix0, ix1], [ix0, ix1], color='gray', linestyle='--', linewidth=0.7, zorder=1)
            x0, x1, y0, y1 = main_lims
            axI.add_patch(mpl.patches.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor='black',
                                                linewidth=0.4, linestyle=':', zorder=2))
            axI.set_xlim(ix0, ix1); axI.set_ylim(iy0, iy1)
            axI.set_xticks(ixt); axI.set_yticks(iyt)
            axI.tick_params(labelsize=5.5, length=2, width=0.4, pad=1)
            for sp in axI.spines.values():
                sp.set_linewidth(0.5)
            axI.set_facecolor('white')
            axI.set_zorder(10)
    # legend: one boxed row above all three panels (inside any panel it would
    # cover points); eta' colours, then plate configuration and overturned
    eta_h = [Line2D([], [], marker='o', color=colors[e], markeredgecolor='black', markeredgewidth=0.25,
                    linestyle='', markersize=4, label="η' = %d" % e) for e in ETAS]
    bc_h = [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
                   linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
           [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray',
                   linestyle='', markersize=4, label='overturned')]
    if Hsize_fig:
        bc_h += [Line2D([], [], marker='o', color='gray', markeredgecolor='black', markeredgewidth=0.25,
                        linestyle='', markersize=np.sqrt(MEAN_S[h]) * 0.75,
                        label=('H = %d km' if h == HS[0] else '%d km') % h) for h in HS]
    plt.tight_layout()
    if strip_fig:
        nH = len(HS) if Hsize_fig else 0
        for h_, e in zip(eta_h, ETAS):
            h_.set_label('%d' % e)
        for h_, h in zip(bc_h[len(bc_h) - nH:], HS):
            h_.set_label('%d km' % h)
        header = lambda t: Line2D([], [], linestyle='', marker='', label=t)
        bc_h[1].set_label('free')          # short labels, so the strip fits the panels' width
        hs = [header("viscosity ratio η':")] + eta_h + [header('plates:')] + bc_h[:len(bc_h) - nH]
        if Hsize_fig:
            hs += [header('thickness:')] + bc_h[len(bc_h) - nH:]
        leg = fig.legend(handles=hs, ncol=len(hs), loc='lower center', bbox_to_anchor=(0.5, 0.91),
                         fontsize=6.5, frameon=False, borderpad=0.5, handletextpad=0.1, handlelength=1.0,
                         columnspacing=0.5)
        for h_, t in zip(hs, leg.get_texts()):
            if h_.get_marker() in ('', 'None', None):
                t.set_fontweight('bold')
        # centre the strip over the panel frames, just above the top tick labels
        fig.canvas.draw()
        r_ = fig.canvas.get_renderer()
        inv = fig.transFigure.inverted()
        x_c = 0.5 * (axP.get_position().x0 + axL.get_position().x1)
        top = max(inv.transform(a.get_tightbbox(r_))[1, 1] for a in (axP, axQ, axL))
        leg.set_bbox_to_anchor((x_c, top + 0.03), transform=fig.transFigure)
        fig.canvas.draw()
        # the box is drawn here rather than as the legend frame: the frame
        # counts the empty handle slot before the first header, which puts the
        # content right of centre.  Box = visible content + equal padding,
        # content shifted so the box is centred over the panels
        lb = leg.get_window_extent()
        tx = [t.get_window_extent(r_) for t in leg.get_texts()]
        cx0, cx1 = tx[0].x0, tx[-1].x1
        shift = (x_c - inv.transform(((cx0 + cx1) / 2, 0))[0])
        leg.set_bbox_to_anchor((x_c + shift, top + 0.03), transform=fig.transFigure)
        fig.canvas.draw()
        tx = [t.get_window_extent(r_) for t in leg.get_texts()]
        hb = [h.get_window_extent(r_) for h in leg.legend_handles]
        cx0, cx1 = tx[0].x0, tx[-1].x1
        cy0 = min(b.y0 for b in tx + hb if b.height > 0); cy1 = max(b.y1 for b in tx + hb if b.height > 0)
        pad = 0.5 * 6.5 / 72 * fig.dpi       # borderpad of 0.5 font sizes
        (bx0, by0), (bx1, by1) = inv.transform(((cx0 - pad, cy0 - pad), (cx1 + pad, cy1 + pad)))
        fig.patches.append(mpl.patches.Rectangle((bx0, by0), bx1 - bx0, by1 - by0, transform=fig.transFigure,
                                                 facecolor='white', edgecolor='black', linewidth=0.5, zorder=-1))
        lb = mpl.transforms.Bbox.from_extents(cx0 - pad, cy0 - pad, cx1 + pad, cy1 + pad)
        span = (axP.get_tightbbox(r_).x0, axL.get_tightbbox(r_).x1)
        print('strip legend %.2f in wide, panels (with labels) %.2f in' % (lb.width / fig.dpi, (span[1] - span[0]) / fig.dpi))
        plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
        plt.savefig(plot_name + '.pdf', bbox_inches='tight')
        print('wrote', plot_name + '.png/.pdf')
        sys.exit(0)
    leg_y = 0.91
    leg = fig.legend(handles=eta_h + bc_h, ncol=len(eta_h + bc_h), loc='lower center',
                     bbox_to_anchor=(0.5, leg_y), fontsize=6.5, frameon=True, fancybox=False,
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
