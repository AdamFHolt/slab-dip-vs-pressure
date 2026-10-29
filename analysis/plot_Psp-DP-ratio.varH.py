#!/bin/python
"""
Supplement figure S10, redrawn for the variable-thickness suite: the share of the
across-slab pressure difference carried by the sub-slab side, P_subslab / DP,
split by plate configuration.

Same selection as Fig. 5 (plot_forces.varH.py ... pooled nodrho71 alleta): every
H = 70/100/130 model, eta' = 50-1000, outputs from t = 11, coverage flag, the
fixed-SP 50 trim, flipped midplanes (dip > 180) dropped, overturned models out.

Reads text_files/withT_varH/, whose column 15 is P_subslab already referenced to
the far-field pressure at the left wall at the analysis depth (extract_properties.py
subtracts get_farfieldP), so no text_files/Pref/ file is needed.  The submitted
S10 (plot_Psp-DP-ratio.py) read the older withT_covered/ files, which hold the raw
pressure, and subtracted the Pref/ file -- but for the five fixed-SP models it
subtracted the fixed-OP model's Pref, which put their ratios at ~0.4 instead of
~0.2 (checked 2026-10-29: each model's own reference reproduces the values here).

The ratio depends on where the reference pressure is taken (here the left wall,
beneath the subducting plate, at the analysis depth); DP itself does not.

Usage:  python3 plot_Psp-DP-ratio.varH.py 300.e3 10.e3 10.e3 1.e3
"""
import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from functions import nominal_plate_thickness

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams['font.size'] = 7
mpl.rcParams['axes.labelsize'] = 7.5
mpl.rcParams['xtick.labelsize'] = 6.5
mpl.rcParams['ytick.labelsize'] = 6.5

analysis_depth    = float(sys.argv[1])
analysis_depth_dz = float(sys.argv[2])
ds                = float(sys.argv[3])
dz                = float(sys.argv[4])

tmin = 11 - 8
suffix = '.z%s.shear-dz%s.ds%s.prof-dz%skm' % (analysis_depth/1.e3, analysis_depth_dz/1.e3, ds/1.e3, dz/1.e3)
out_dir = 'plots/DP-comparisons/compilations/manuscript_figs/'
os.makedirs(out_dir, exist_ok=True)
plot_name = out_dir + 'Psp-DP-ratio_supp' + suffix + '.tmin%d' % tmin

DP_ind = 3; Psp_ind = 15; cov_ind = 38

# model names, overturned set and trims as in plot_forces.varH.py
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
    if H == 100:
        return H100[bc][eta]
    return '2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % ('' if bc == 'free' else bc + '_', eta, H)

BCS  = ['FixedSP', 'free', 'FixedOP']
ETAS = [50, 250, 375, 500, 1000]
HS   = [70, 100, 130]
OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130),
              ('FixedOP', 500, 100), ('FixedOP', 1000, 100),
              ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70),
              ('FixedOP', 375, 100)}
LAST_OUTPUT = {('FixedSP', 50, 100): 22}
bc_label = {'free': 'free plates', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
bc_color = {'FixedSP': '#3b6e8f', 'free': '#9a9a9a', 'FixedOP': '#b5651d'}

def load(bc, eta, H):
    f = 'text_files/withT_varH/' + name(bc, eta, H) + suffix + '.txt'
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f)[tmin:]
    d = d[d[:, cov_ind] == 1]
    d = d[d[:, 0] <= LAST_OUTPUT.get((bc, eta, H), 1e9)]
    d = d[d[:, 5] <= 180]
    r = d[:, Psp_ind] / d[:, DP_ind]
    return r[np.isfinite(r)]

ratios = {bc: [] for bc in BCS}
means  = {bc: [] for bc in BCS}
print('%-8s %-5s %-4s %4s %7s' % ('plates', 'eta', 'H', 'N', 'median'))
for bc in BCS:
    for e in ETAS:
        for H in HS:
            if (bc, e, H) in OVERTURNED:
                continue
            r = load(bc, e, H)
            if r is None or len(r) == 0:
                continue
            ratios[bc].append(r); means[bc].append(np.median(r))
            print('%-8s %-5d %-4d %4d %7.2f' % (bc, e, H, len(r), np.median(r)))
ratios = {bc: np.concatenate(v) for bc, v in ratios.items()}
allr = np.concatenate(list(ratios.values()))
print()
for bc in BCS + ['all']:
    r = allr if bc == 'all' else ratios[bc]
    print('%-8s N %3d  p10 %.2f  median %.2f  p90 %.2f  outside -0.1..1.1: %d'
          % (bc, len(r), *np.percentile(r, [10, 50, 90]), np.sum((r < -0.1) | (r > 1.1))))
for H in HS:
    r = np.concatenate([load(bc, e, H) for bc in BCS for e in ETAS
                        if (bc, e, H) not in OVERTURNED and load(bc, e, H) is not None])
    print('H = %d    N %3d  median %.2f' % (H, len(r), np.median(r)))

fig, ax = plt.subplots(figsize=(3.4, 2.4))
bins = np.arange(-0.1, 1.1 + 1e-9, 0.05)
for bc in BCS:
    r = np.clip(ratios[bc], -0.1, 1.1)
    ax.hist(r, bins=bins, color=bc_color[bc], alpha=0.35, zorder=2)
    ax.hist(r, bins=bins, histtype='step', color=bc_color[bc], lw=0.9, zorder=3,
            label='%s (median %.2f)' % (bc_label[bc], np.median(ratios[bc])))
    ax.axvline(np.median(ratios[bc]), color=bc_color[bc], lw=0.9, ls='--', zorder=1.5)
ax.axvline(0.5, color='lightgray', lw=1, zorder=0)
ax.set_xlim(-0.1, 1.1)
ax.set_xlabel(r'$P_{subslab}\, /\, \Delta P$')
ax.set_ylabel('number of model outputs')
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
ax.tick_params(which='both', direction='out')
leg = ax.legend(loc='upper right', fontsize=6.5, frameon=True, fancybox=False, edgecolor='black', framealpha=1,
                handlelength=1.5)
leg.get_frame().set_linewidth(0.5)
fig.savefig(plot_name + '.png', bbox_inches='tight', dpi=600)
fig.savefig(plot_name + '.pdf', bbox_inches='tight')
print('wrote', plot_name + '.png/pdf')
