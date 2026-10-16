#!/bin/python
"""
Supplement figure: how the bending length L_K of the n = 3 scaling is fitted.

dQ_s/ds = eta*H^3*K*v_c / (3*L_K^3) is fitted through the origin (least squares
in MPa) to every normal output of the Fig. 5 selection: text_files/withT_varH/,
t >= 11, coverage flag, the fixed-SP eta' = 50 trim, all five eta' at H = 100,
no drho71 run, overturned models left out.  The models, the overturned set and
the trim are read from plot_forces.varH.py so the two cannot drift apart.

  (a) rms misfit against L_K at 250/300/350 km, H = the nominal plate thickness
      (as in Fig. 5 and the Earth maps)
  (b) the best-fitting L_K for all models and for subsets (H, plate
      configuration) at the three depths, with 16-84 % ranges from a bootstrap
      over models; the H = 70 and fixed-OP subsets, whose dQ/ds is ~0, give no
      positive fit and are not drawn

Writes plots/DP-comparisons/compilations/manuscript_figs/LK_fit.tmin3.png/pdf
and prints the numbers for the caption.

Usage:  python3 plot_LK_fit.py [L_adopted_km]      (default 200)
"""
import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from functions import nominal_plate_thickness

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})

L_adopted = float(sys.argv[1]) if len(sys.argv) > 1 else 200.

# model names, OVERTURNED and LAST_OUTPUT from plot_forces.varH.py (that script
# runs on import, so only these blocks are executed)
src = open('plot_forces.varH.py').read()
ns = {}
for blk in ('H100 = {', 'def name(', 'OVERTURNED = {', 'LAST_OUTPUT = {'):
    i = src.index(blk)
    exec(src[i:src.index('\n\n', i)], ns)

K_ind, sn_ind, vc_ind, cov_ind = 11, 6, 19, 38
mant_visc  = 2.5e20
cmyr_to_ms = 0.01/(365.25*24*3600)
tmin = 3                                   # rows from t = 11
DEPTHS = (250, 300, 350)
BCS, ETAS, HS = ('free', 'FixedSP', 'FixedOP'), (50, 250, 375, 500, 1000), (70, 100, 130)

def load(z):
    """per normal model: X = eta*K*v_c*H^3 (nominal H) and dQ/ds, SI"""
    out = {}
    for bc in BCS:
        for e in ETAS:
            for H in HS:
                if (bc, e, H) in ns['OVERTURNED']:
                    continue
                m = ns['name'](bc, e, H)
                f = 'text_files/withT_varH/%s.z%.1f.shear-dz10.0.ds10.0.prof-dz1.0km.txt' % (m, z)
                if not os.path.exists(f):
                    continue
                d = np.loadtxt(f)[tmin:]
                d = d[d[:, cov_ind] == 1]
                d = d[d[:, 0] <= ns['LAST_OUTPUT'].get((bc, e, H), 1e9)]
                out[(bc, e, H)] = dict(Xn=e * mant_visc * d[:, K_ind] * d[:, vc_ind] * cmyr_to_ms
                                          * (nominal_plate_thickness(m)*1e3)**3, y=-d[:, sn_ind])
    return out

def c_fit(models, key):
    x = np.concatenate([d[key] for d in models]); y = np.concatenate([d['y'] for d in models])
    return np.sum(x*y) / np.sum(x*x)

def L_of(c):                               # km; c = 1/(3 L^3)
    return (1/(3*c))**(1/3) / 1e3 if c > 0 else np.nan

def rms(models, key, L_km):
    x = np.concatenate([d[key] for d in models]); y = np.concatenate([d['y'] for d in models])
    return np.sqrt(np.mean((y - x/(3*(L_km*1e3)**3))**2)) / 1e6

def boot(models, key, nb=2000, seed=0):
    """16-84 % range of L over resamplings of whole models (c <= 0 -> L unbounded)"""
    rng = np.random.default_rng(seed)
    cs = [c_fit([models[i] for i in rng.integers(0, len(models), len(models))], key) for _ in range(nb)]
    c16, c84 = np.percentile(cs, [16, 84])
    return L_of(c84), (L_of(c16) if c16 > 0 else np.inf)

SUBSETS = [('all', lambda k: True), ('H = 100', lambda k: k[2] == 100), ('H = 130', lambda k: k[2] == 130),
           ('free', lambda k: k[0] == 'free'), ('fixed SP', lambda k: k[0] == 'FixedSP'),
           ('H = 70', lambda k: k[2] == 70), ('fixed OP', lambda k: k[0] == 'FixedOP')]
data = {z: load(z) for z in DEPTHS}
Ls = np.arange(120., 351., 1.)
dcol = {250: '#6baed6', 300: '#08519c', 350: '#e6550d'}
dmrk = {250: 'v', 300: 'o', 350: '^'}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw={'width_ratios': [1, 1.25]})

print('%-5s %-9s %4s %3s %6s  %-11s %8s %8s' % ('z', 'subset', 'N', 'mod', 'L_K', '16-84%', 'rms@fit', 'rms@%d' % L_adopted))
fits = {}
for z in DEPTHS:
    for lab, sel in SUBSETS:
        models = [d for k, d in data[z].items() if sel(k)]
        N = sum(len(d['y']) for d in models)
        L = L_of(c_fit(models, 'Xn')); lo, hi = boot(models, 'Xn')
        fits[(z, lab)] = (L, lo, hi)
        print('%-5d %-9s %4d %3d %6.0f  %4.0f-%-6.0f %8.2f %8.2f'
              % (z, lab, N, len(models), L, lo, hi, rms(models, 'Xn', L) if L == L else np.nan, rms(models, 'Xn', L_adopted)))

# (a) misfit curves
for z in DEPTHS:
    models = list(data[z].values())
    r = np.array([rms(models, 'Xn', L) for L in Ls])
    axA.plot(Ls, r, color=dcol[z], linewidth=1.2 if z == 300 else 0.9, label='%d km' % z, zorder=3)
    L = fits[(z, 'all')][0]
    axA.plot(L, rms(models, 'Xn', L), dmrk[z], color=dcol[z], markersize=3.5, markeredgecolor='black',
             markeredgewidth=0.3, zorder=4)
    if z == 300:
        lo, hi = fits[(z, 'all')][1:]
        axA.axvspan(lo, hi, color=dcol[z], alpha=0.12, linewidth=0, zorder=1)
axA.axvline(L_adopted, color='gray', linestyle=':', linewidth=1, zorder=2)
axA.set_xlim(Ls[0], Ls[-1]); axA.set_ylim(2, 6)
axA.set_xlabel(r'$L_K$   [km]'); axA.set_ylabel(r'rms misfit of $-dQ_s/ds$   [MPa]')
axA.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
leg = axA.legend(loc='upper right', fontsize=6.5, frameon=True, fancybox=False, edgecolor='black',
                 facecolor='white', framealpha=1, handlelength=1.8, borderpad=0.5)
leg.get_frame().set_linewidth(0.5)

# (b) L_K by subset and depth
shown = [s for s, _ in SUBSETS if any(fits[(z, s)][0] == fits[(z, s)][0] for z in DEPTHS)
         and s not in ('H = 70', 'fixed OP')]
for i, lab in enumerate(shown):
    for j, z in enumerate(DEPTHS):
        L, lo, hi = fits[(z, lab)]
        if not L == L:
            continue
        x = i + (j - 1) * 0.22
        hi_d = min(hi, 345)
        axB.plot([x, x], [lo, hi_d], color=dcol[z], linewidth=1.2, zorder=2)
        if hi > 345:
            axB.annotate('', xy=(x, 350), xytext=(x, 330), arrowprops=dict(arrowstyle='->', color=dcol[z], lw=1))
        axB.plot(x, L, dmrk[z], color=dcol[z], markersize=4, markeredgecolor='black', markeredgewidth=0.3,
                 zorder=3, label='%d km' % z if i == 0 else None)
axB.axhline(L_adopted, color='gray', linestyle=':', linewidth=1, zorder=1)
axB.set_xticks(range(len(shown)))
axB.set_xticklabels(shown)
axB.set_xlim(-0.5, len(shown) - 0.5); axB.set_ylim(150, 350)
axB.set_ylabel(r'best-fitting $L_K$   [km]')
axB.grid(True, axis='y', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
leg = axB.legend(loc='upper right', fontsize=6.5, ncol=3, frameon=True, fancybox=False, edgecolor='black',
                 facecolor='white', framealpha=1, handletextpad=0.1, columnspacing=0.8, borderpad=0.5)
leg.get_frame().set_linewidth(0.5)

for ax, lab in ((axA, 'a'), (axB, 'b')):
    ax.text(0.04, 0.96, lab + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10,
            bbox=dict(facecolor='white', edgecolor='none', pad=1))

plt.tight_layout()
out_dir = 'plots/DP-comparisons/compilations/manuscript_figs/'
os.makedirs(out_dir, exist_ok=True)
plot_name = out_dir + 'LK_fit.tmin%d' % tmin
plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
print('wrote', plot_name + '.png/.pdf')
