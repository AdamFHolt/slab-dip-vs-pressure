"""
Summary of Test A (extract_testA_kinematics.py): per-model distribution of
r = -(dv_s/dn)/omega_t, the stress-based ratio tau_sn/(eta omega_t), and the
extraction check tau_sn/(eta 2eps_sn).  Also the share of omega_t carried by
the steady translation term -K v_s, to separate "translating through a fixed
bend" from a rotating slab.

Figure: plots/tmp/testA_kinematics.z<zkm>.png (new file).
Usage:  python3 plot_testA_kinematics.py [analysis_depth_m=300e3]
"""
import glob
import os
import re
import sys

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
            b = b[len(p):]
            break
    return re.sub(r'\.z[\d.]+\.txt$', '', b)


def eta_of(lab):
    m = re.search(r'(\d+)plates', lab)
    return int(m.group(1)) if m else 500


def H_of(lab):
    m = re.search(r'_H(\d+)km', lab)
    return int(m.group(1)) if m else 100


def bc_of(lab):
    return 'fixedSP' if 'FixedSP' in lab else 'fixedOP' if 'FixedOP' in lab else 'free'


files = sorted(glob.glob('text_files/testA/*.z%s.txt' % zkm))
rows = []
print('%-28s %4s %6s %13s %12s %12s %10s' % ('model', 'N', 'r_med', 'r_IQR', 'tau/(eta w)', 'tau/(eta2e)', '-Kvs/w'))
allr = {}
for f in files:
    d = np.loadtxt(f)
    if d.ndim == 1:
        d = d[None, :]
    d = d[d[:, 26] == 1]                       # coverage flag
    lab = label(f)
    r = d[:, 10]
    kvs = -d[:, 2] * d[:, 6] * 0.01 / (365.25 * 24 * 3600) / d[:, 9]   # (-K v_s)/omega_t
    chk = d[:, 16]
    chk = chk[np.abs(d[:, 12]) > 0.2 * np.abs(d[:, 9])]   # skip rows where 2eps ~ 0
    tag = ' OT' if lab in OVERTURNED else ''
    print('%-28s %4d %6.2f %6.2f-%-6.2f %12.2f %12.2f %10.2f%s'
          % (lab, len(r), np.median(r), np.percentile(r, 25), np.percentile(r, 75),
             np.median(d[:, 15]), np.median(np.abs(chk)) if len(chk) else np.nan, np.median(kvs), tag))
    allr[lab] = (d[:, 0], r, d[:, 15], eta_of(lab), H_of(lab), bc_of(lab), lab in OVERTURNED)
    rows.append(d)

allrows = np.vstack(rows)
normal = np.concatenate([v[1] for k, v in allr.items() if not v[6]])
print('\nall normal models: N=%d  r median %.2f, IQR %.2f-%.2f, 10-90%% %.2f-%.2f; fraction r>0.7: %.2f, r<0.3: %.2f'
      % (len(normal), np.median(normal), *np.percentile(normal, [25, 75]), *np.percentile(normal, [10, 90]),
         np.mean(normal > 0.7), np.mean(normal < 0.3)))
for H in (70, 100, 130):
    rr = np.concatenate([v[1] for v in allr.values() if v[4] == H and not v[6]])
    print('  H=%3d: N=%3d  r median %.2f  IQR %.2f-%.2f' % (H, len(rr), np.median(rr), *np.percentile(rr, [25, 75])))

# --- figure
colors = {50: 'tan', 250: 'peru', 375: 'firebrick', 500: 'maroon', 1000: 'black'}
markers = {'free': 'o', 'fixedSP': 'v', 'fixedOP': '^'}
sizes = {70: 12, 100: 20, 130: 30}
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
for lab, (t, r, ratio, e, H, bc, ot) in allr.items():
    kw = dict(color=colors[e], marker=markers[bc], s=sizes[H], linewidth=0.5)
    if ot:
        kw.update(facecolors='none', edgecolor=colors[e], alpha=0.6)
    else:
        kw.update(edgecolor='black')
    axs[0].scatter(t, r, **kw)
    axs[2].scatter(r, ratio, **kw)
for ax in axs[:1]:
    ax.axhline(1, color='gray', ls='--', lw=1); ax.axhline(0, color='gray', ls='--', lw=1)
axs[0].set_xlabel('output number'); axs[0].set_ylabel(r'$r = -(\partial v_s/\partial n)\,/\,\omega_t$')
axs[0].set_ylim(-0.5, 2)
axs[0].set_title('r = 1: Kirchhoff (bending);  r = 0: shear mode (Eq. 6)', fontsize=9)
axs[1].hist(normal, bins=np.linspace(-0.5, 2, 51), color='maroon', alpha=0.8)
axs[1].axvline(1, color='gray', ls='--'); axs[1].axvline(0, color='gray', ls='--')
axs[1].set_xlabel('r (normal models, all timesteps)'); axs[1].set_ylabel('count')
axs[2].axhline(1, color='gray', ls='--', lw=1); axs[2].axhline(0, color='gray', ls='--', lw=1)
axs[2].set_xlim(-0.5, 2); axs[2].set_ylim(-0.5, 1.5)
axs[2].set_xlabel('r'); axs[2].set_ylabel(r'$\bar\tau_{sn}\,/\,(\eta\,\omega_t)$  (Test A$^\prime$)')
from matplotlib.lines import Line2D
h = [Line2D([], [], marker='o', color=colors[e], ls='', label=r"$\eta'$=%d" % e) for e in (50, 250, 375, 500, 1000)]
h += [Line2D([], [], marker=markers[b], color='gray', ls='', label=b) for b in markers]
h += [Line2D([], [], marker='o', color='gray', ls='', markersize=np.sqrt(sizes[H]), label='H=%d' % H) for H in sizes]
axs[0].legend(handles=h, fontsize=6, ncol=3, loc='lower right', frameon=False)
plt.tight_layout()
out = 'plots/tmp/testA_kinematics.z%s.png' % zkm
plt.savefig(out, dpi=150)
print('wrote', out)
