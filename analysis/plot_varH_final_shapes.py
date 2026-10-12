#!/usr/bin/env python3
"""
Slab shapes at the last analysed output of every model in the variable-H
compilation (free, fixed-SP, fixed-OP x H = 70/100/130 x eta' = 250/500/1000,
plus the B-matched drho71 run), to classify them as normal or overturned.

Each panel overlays the c_llith = 0.5 contour (5 km grid, x relative to where
the slab top crosses 300 km depth) for the three eta' at one (BC, H).  The
printed table gives, for the deepest 600 km of the slab, the horizontal offset
of the tip from the slab's furthest-right point: a tip well to the LEFT of that
point (negative offset) means the slab folded back toward the SP, i.e. it is
overturned.  Classification remains a judgement call from the figure.

Usage: python3 plot_varH_final_shapes.py   (from analysis/)
Writes plots/tmp/varH_final_shapes.png
"""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from multiprocessing import Pool

P = '2D_compositional_subd_lower-res_new_'
# (bc, eta, H): (model, maxt)   -- last analysed output is maxt-1
MODELS = {
 ('free',250,70): (P+'250plates_H70km', 39), ('free',500,70): (P+'500plates_H70km', 32), ('free',1000,70): (P+'1000plates_H70km', 31),
 ('free',250,100): (P+'250plates', 28), ('free',500,100): ('2D_compositional_subd_lower-res_new2', 31), ('free',1000,100): (P+'1000plates', 30),
 ('free',250,130): (P+'250plates_H130km', 26), ('free',500,130): (P+'500plates_H130km', 30), ('free',1000,130): (P+'1000plates_H130km', 33),
 ('FixedSP',250,70): (P+'FixedSP_250plates_H70km', 34), ('FixedSP',500,70): (P+'FixedSP_500plates_H70km', 32), ('FixedSP',1000,70): (P+'FixedSP_1000plates_H70km', 30),
 ('FixedSP',250,100): (P+'FixedSP_250plates', 30), ('FixedSP',500,100): ('2D_compositional_subd_FixedSP_lower-res_new2', 30), ('FixedSP',1000,100): (P+'FixedSP_1000plates2', 30),
 ('FixedSP',250,130): (P+'FixedSP_250plates_H130km', 46), ('FixedSP',500,130): (P+'FixedSP_500plates_H130km', 42),
 ('FixedOP',250,70): (P+'FixedOP_250plates_H70km', 36), ('FixedOP',500,70): (P+'FixedOP_500plates_H70km', 33), ('FixedOP',1000,70): (P+'FixedOP_1000plates_H70km', 36),
 ('FixedOP',250,100): (P+'FixedOP_250plates', 28), ('FixedOP',500,100): ('2D_compositional_subd_FixedOP_lower-res_new', 28), ('FixedOP',1000,100): (P+'FixedOP_1000plates', 30),
 ('FixedOP',250,130): (P+'FixedOP_250plates_H130km', 33), ('FixedOP',500,130): (P+'FixedOP_500plates_H130km', 40), ('FixedOP',1000,130): (P+'FixedOP_1000plates_H130km', 42),
 ('drho71',500,70): (P+'500plates_H70km_drho71', 39),
}
YMAX = 1450e3

def contour(key):
    m, maxt = MODELS[key]; t = maxt - 1
    fn = 'csv_outputs/%s/full.%d.csv' % (m, t)
    hdr = open(fn).readline().strip().replace('"', '').split(',')
    d = np.loadtxt(fn, delimiter=',', skiprows=1)
    xc, yc, lc = hdr.index('Points:0'), hdr.index('Points:1'), hdr.index('llith')
    X, Y = np.meshgrid(np.arange(1000e3, 5000e3 + 1, 5e3), np.arange(0, YMAX + 1, 5e3))
    ll = griddata((d[:, xc], d[:, yc]), d[:, lc], (X, Y), method='linear')
    cs = plt.contour(X / 1e3, (YMAX - Y) / 1e3, ll, levels=[0.5])
    seg = max(cs.allsegs[0], key=len)          # the SP + slab contour
    plt.close('all')
    slab = seg[seg[:, 1] > 150]                 # below the plates
    x_ref = seg[np.abs(seg[:, 1] - 300) < 5, 0].max()   # slab top at 300 km (the plates end above 150 km)
    tip = slab[np.argmax(slab[:, 1])]
    deep = slab[slab[:, 1] > tip[1] - 600]
    offset = tip[0] - deep[:, 0].max()          # km; strongly negative = folded back
    return key, seg, x_ref, tip[1], offset

if __name__ == '__main__':
    with Pool(14) as pool:
        res = pool.map(contour, list(MODELS))
    out = {r[0]: r for r in res}
    bcs = ['free', 'FixedSP', 'FixedOP']; Hs = [70, 100, 130]
    col = {250: 'tab:blue', 500: 'tab:orange', 1000: 'tab:red'}
    fig, axs = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
    print('%-8s %5s %4s  %-45s %7s %10s' % ('bc', 'eta', 'H', 'model (last output)', 'tip km', 'tip offset km'))
    for i, bc in enumerate(bcs):
        for j, H in enumerate(Hs):
            ax = axs[i, j]
            for eta in [250, 500, 1000]:
                key = (bc, eta, H)
                if key not in out: continue
                _, seg, x_ref, ztip, off = out[key]
                ax.plot(seg[:, 0] - x_ref, seg[:, 1], color=col[eta], lw=1, label=r"$\eta'$=%d" % eta)
                print('%-8s %5d %4d  %-45s %7.0f %10.0f' % (bc, eta, H, MODELS[key][0].replace('2D_compositional_subd_', '') + ' (%d)' % (MODELS[key][1] - 1), ztip, off))
            if bc == 'free' and H == 70:
                _, seg, x_ref, ztip, off = out[('drho71', 500, 70)]
                ax.plot(seg[:, 0] - x_ref, seg[:, 1], color='tab:orange', lw=1, ls='--', label='500, drho71')
                print('%-8s %5d %4d  %-45s %7.0f %10.0f' % ('drho71', 500, 70, MODELS[('drho71', 500, 70)][0].replace('2D_compositional_subd_', '') + ' (38)', ztip, off))
            ax.axhline(660, color='gray', lw=0.5, ls=':')
            ax.set_xlim(-800, 1200); ax.set_ylim(1450, 0); ax.set_aspect('equal')
            ax.set_title('%s, H = %d km' % (bc, H), fontsize=10)
            if i == 2: ax.set_xlabel('x - x(slab top at 300 km)  [km]')
            if j == 0: ax.set_ylabel('depth [km]')
            ax.legend(fontsize=7, loc='lower right')
    plt.tight_layout(); plt.savefig('plots/tmp/varH_final_shapes.png', dpi=150)
    print('wrote plots/tmp/varH_final_shapes.png')
