"""
Summary of Test C (extract_testC_moment.py), per model and pooled:

  identity   (dM/ds + C) / Q_s            should be ~1 in any regime
  partition  dM/ds / Q_s  and  C / Q_s    bending: moment gradient carries Q_s;
                                          shear mode: the mantle couple does
  dQ/ds      fitted vs extraction (col 13), and d2M/ds2 + dC/ds vs dQ/ds
  moment     M / (eta H^3 kdot/3)         thin-sheet prefactor check
  kdot       kinematic / stress           should be ~1

Figure: plots/tmp/testC_moment.z<zkm>.png
Usage:  python3 plot_testC_moment.py [analysis_depth_m=300e3]
"""
import glob, os, re, sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

depth = float(sys.argv[1]) if len(sys.argv) > 1 else 300.e3
zkm = '%.1f' % (depth / 1.e3)
OVERTURNED = {'new_1000plates', 'new_FixedOP_1000plates', 'FixedOP_lower-res_new',
              'new_FixedOP_375plates', 'new_500plates_H130km', 'new_1000plates_H130km'}


def label(f):
    b = os.path.basename(f)
    for p in ('2D_compositional_subd_lower-res_', '2D_compositional_subd_'):
        if b.startswith(p):
            b = b[len(p):]; break
    return re.sub(r'\.z[\d.]+\.txt$', '', b)


def H_of(lab):
    m = re.search(r'_H(\d+)km', lab); return int(m.group(1)) if m else 100


def med(x):
    x = x[np.isfinite(x)]
    return np.median(x) if len(x) else np.nan


files = sorted(glob.glob('text_files/testC/*.z%s.txt' % zkm))
print('%-28s %3s | %8s %8s %6s | %7s %7s %7s | %6s %6s' % (
    'model', 'N', 'ident', 'dM/ds/Q', 'C/Q', 'dQfit', 'dQtxt', 'd2M+dC', 'M/Mest', 'kd k/s'))
pool = {}
for f in files:
    d = np.loadtxt(f)
    if d.ndim == 1:
        d = d[None, :]
    d = d[d[:, 6] == 1]
    lab = label(f)
    # use |Q_s| > 100 GN/m rows for the ratios, so tiny Q does not blow them up
    big = np.abs(d[:, 7]) > 1e11
    ident = d[big, 11] / d[big, 7]
    part_M = d[big, 9] / d[big, 7]
    part_C = d[big, 10] / d[big, 7]
    tag = ' OT' if lab in OVERTURNED else ''
    print('%-28s %3d | %8.2f %8.2f %6.2f | %7.2f %7.2f %7.2f | %6.2f %6.2f%s' % (
        lab, len(d), med(ident), med(part_M), med(part_C),
        med(d[:, 12]) / 1e6, med(d[:, 13]) / 1e6, med(d[:, 16]) / 1e6,
        med(d[:, 17] / d[:, 23]), med(d[:, 21] / d[:, 22]), tag))
    pool[lab] = (d, H_of(lab), lab in OVERTURNED, big)

normal = np.vstack([v[0] for v in pool.values() if not v[2]])
bigN = np.abs(normal[:, 7]) > 1e11
print('\nnormal models pooled (N=%d, %d with |Q_s|>100 GN/m):' % (len(normal), bigN.sum()))
print('  identity (dM/ds+C)/Q: median %.2f, IQR %.2f-%.2f' % (med(normal[bigN, 11] / normal[bigN, 7]),
      *np.nanpercentile(normal[bigN, 11] / normal[bigN, 7], [25, 75])))
print('  partition dM/ds / Q: median %.2f ;  C / Q: median %.2f' % (
      med(normal[bigN, 9] / normal[bigN, 7]), med(normal[bigN, 10] / normal[bigN, 7])))
ok = np.isfinite(normal[:, 16]) & np.isfinite(normal[:, 13])
print('  r(dQ/ds txt, dQ/ds fit) = %.2f ; r(dQ/ds txt, d2M/ds2 + dC/ds) = %.2f ; r(dQ/ds txt, d2M/ds2) = %.2f'
      % (np.corrcoef(normal[ok, 13], normal[ok, 12])[0, 1], np.corrcoef(normal[ok, 13], normal[ok, 16])[0, 1],
         np.corrcoef(normal[ok, 13], normal[ok, 14])[0, 1]))
print('  M / (eta H^3 kdot/3): median %.2f ;  kdot kin/stress: median %.2f' % (
      med(normal[:, 17] / normal[:, 23]), med(normal[:, 21] / normal[:, 22])))
for H in (70, 100, 130):
    sel = np.array([H_of(l) == H for l in []])  # placeholder
    rows = [v[0] for v in pool.values() if v[1] == H and not v[2]]
    if not rows:
        continue
    dH = np.vstack(rows); b = np.abs(dH[:, 7]) > 1e11
    print('  H=%3d: N=%3d  med Q_s %6.0f GN/m  med dM/ds %6.0f  med C %6.0f  | med dQ/ds txt %5.2f  fit %5.2f  d2M+dC %5.2f MPa | med |M| %.2e N'
          % (H, len(dH), med(dH[:, 7]) / 1e9, med(dH[:, 9]) / 1e9, med(dH[:, 10]) / 1e9,
             med(dH[:, 13]) / 1e6, med(dH[:, 12]) / 1e6, med(dH[:, 16]) / 1e6, med(np.abs(dH[:, 17]))))

# figure
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
cols = {70: 'tab:blue', 100: 'tab:gray', 130: 'tab:red'}
for lab, (d, H, ot, big) in pool.items():
    kw = dict(s=14, color=cols[H], edgecolor='k', linewidth=0.3, zorder=3)
    if ot:
        kw.update(facecolors='none', edgecolor=cols[H], alpha=0.6)
    axs[0].scatter(d[:, 7] / 1e9, d[:, 11] / 1e9, **kw)
    axs[1].scatter(d[:, 7] / 1e9, d[:, 9] / 1e9, **kw)
    axs[1].scatter(d[:, 7] / 1e9, d[:, 10] / 1e9, marker='x', s=10, color=cols[H], alpha=0.6)
    axs[2].scatter(d[:, 13] / 1e6, d[:, 16] / 1e6, **kw)
for ax, (xl, yl, lim) in zip(axs, [(r'$Q_s$ [GN/m]', r'$dM/ds + C$ [GN/m]', (-500, 2500)),
                                    (r'$Q_s$ [GN/m]', r'$dM/ds$ (dots),  $C$ (crosses) [GN/m]', (-500, 2500)),
                                    (r'measured $dQ_s/ds$ [MPa]', r'$d^2M/ds^2 + dC/ds$ [MPa]', (-10, 30))]):
    ax.plot(lim, lim, 'k--', lw=0.8); ax.axhline(0, color='0.8', lw=0.8); ax.axvline(0, color='0.8', lw=0.8)
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(color='0.92')
from matplotlib.lines import Line2D
axs[0].legend(handles=[Line2D([], [], marker='o', ls='', color=cols[H], label='H = %d' % H) for H in cols],
              fontsize=7, frameon=False)
axs[0].set_title('moment balance identity', fontsize=9)
axs[1].set_title('partition of $Q_s$', fontsize=9)
axs[2].set_title(r'$dQ_s/ds$ from the moment (noisy 2nd derivative)', fontsize=9)
plt.tight_layout()
out = 'plots/tmp/testC_moment.z%s.png' % zkm
plt.savefig(out, dpi=140)
print('wrote', out)
