#!/usr/bin/env python3
"""
Supplement figure: slab shape at the last analysed output of every model in the
variable-H compilation (fixed SP, free, fixed OP x H = 70/100/130 km x
eta' = 250/500/1000), the basis of the normal / overturned classification used
in Figs 3-5 (OVERTURNED in plot_forces.varH.py; overturned runs drawn dashed).
Each panel overlays the c_llith = 0.5 contour (5 km grid, x relative to where
the slab top crosses 300 km depth) for the three eta' at one (plates, H).
The printed table gives, for the deepest 600 km of the slab, the horizontal
offset of the tip from the slab's furthest-right point: a tip well to the LEFT
(negative offset) means the slab folded back toward the SP.  The
classification itself is the judgement call recorded in OVERTURNED.
Usage: python3 plot_varH_final_shapes.py   (from analysis/)
Writes plots/DP-comparisons/compilations/manuscript_figs/final_shapes_supp.png/.pdf
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
 ('FixedSP',250,130): (P+'FixedSP_250plates_H130km', 46), ('FixedSP',500,130): (P+'FixedSP_500plates_H130km', 42), ('FixedSP',1000,130): (P+'FixedSP_1000plates_H130km', 37),
 ('FixedOP',250,70): (P+'FixedOP_250plates_H70km', 36), ('FixedOP',500,70): (P+'FixedOP_500plates_H70km', 33), ('FixedOP',1000,70): (P+'FixedOP_1000plates_H70km', 36),
 ('FixedOP',250,100): (P+'FixedOP_250plates', 28), ('FixedOP',500,100): ('2D_compositional_subd_FixedOP_lower-res_new', 28), ('FixedOP',1000,100): (P+'FixedOP_1000plates', 30),
 ('FixedOP',250,130): (P+'FixedOP_250plates_H130km', 33), ('FixedOP',500,130): (P+'FixedOP_500plates_H130km', 40), ('FixedOP',1000,130): (P+'FixedOP_1000plates_H130km', 42),
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

OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130),
              ('FixedOP', 500, 100), ('FixedOP', 1000, 100),
              ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70)}        # as plot_forces.varH.py (eta' = 375 is not in this compilation)

if __name__ == '__main__':
    import os
    import matplotlib as mpl
    import matplotlib.font_manager as fm
    from matplotlib.lines import Line2D
    font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
    if os.path.exists(font_path):
        fm.fontManager.addfont(font_path)
        mpl.rcParams['font.family'] = 'Myriad Pro'
    mpl.rcParams.update({'font.size': 7, 'axes.labelsize': 9, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})
    with Pool(14) as pool:
        res = pool.map(contour, list(MODELS))
    out = {r[0]: r for r in res}
    bcs = ['FixedSP', 'free', 'FixedOP']; Hs = [70, 100, 130]          # rows as Fig. 2
    bc_label = {'FixedSP': 'fixed SP', 'free': 'free plates', 'FixedOP': 'fixed OP'}
    col = {250: 'peru', 500: 'maroon', 1000: 'black'}                  # the paper's eta' colours
    fig, axs = plt.subplots(3, 3, figsize=(7.0, 5.6), sharex=True, sharey=True)
    print('%-8s %5s %4s  %-45s %7s %10s  %s' % ('bc', 'eta', 'H', 'model (last output)', 'tip km', 'tip offset km', ''))
    n = 0
    for i, bc in enumerate(bcs):
        for j, H in enumerate(Hs):
            ax = axs[i, j]
            for eta in [250, 500, 1000]:
                key = (bc, eta, H)
                if key not in out: continue
                _, seg, x_ref, ztip, off = out[key]
                ot = key in OVERTURNED
                ax.plot(seg[:, 0] - x_ref, seg[:, 1], color=col[eta], lw=0.9, ls=(0, (3, 1.5)) if ot else '-',
                        zorder=3 if ot else 4)
                print('%-8s %5d %4d  %-45s %7.0f %10.0f  %s' % (bc, eta, H, MODELS[key][0].replace('2D_compositional_subd_', '')
                      + ' (%d)' % (MODELS[key][1] - 1), ztip, off, 'overturned' if ot else ''))
            ax.axhline(660, color='gray', lw=0.5, ls=':', zorder=1)
            ax.set_xlim(-800, 1200); ax.set_ylim(1450, 0); ax.set_aspect('equal')
            ax.set_xticks([-500, 0, 500, 1000]); ax.set_yticks([0, 500, 1000])
            ax.tick_params(length=2.5)
            ax.text(0.03, 0.95, 'abcdefghi'[n] + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10,
                    bbox=dict(facecolor='white', edgecolor='none', pad=1))
            ax.text(0.97, 0.95, '%s, H = %d km' % (bc_label[bc], H), transform=ax.transAxes, ha='right', va='top', fontsize=8, zorder=10,
                    bbox=dict(facecolor='white', edgecolor='none', pad=1))
            if i == 0 and j == 0:
                ax.text(-780, 660, '660 km', color='gray', fontsize=6.5, ha='left', va='bottom')
            if i == 2: ax.set_xlabel('$x - x$(slab top at 300 km)   [km]')
            if j == 0: ax.set_ylabel('depth   [km]')
            n += 1
    hs = [Line2D([], [], color=col[e], lw=1.2, label="\u03b7' = %d" % e) for e in [250, 500, 1000]] + \
         [Line2D([], [], color='gray', lw=1.2, ls=(0, (3, 1.5)), label='overturned')]
    leg = axs[0, 0].legend(handles=hs, loc='lower left', fontsize=6.5, frameon=True, fancybox=False, edgecolor='black',
                           facecolor='white', framealpha=1, borderpad=0.5, labelspacing=0.3, handlelength=1.6, borderaxespad=0.5)
    leg.get_frame().set_linewidth(0.5)
    plt.tight_layout(h_pad=0.6, w_pad=0.6)
    out_dir = 'plots/DP-comparisons/compilations/manuscript_figs/'
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(out_dir + 'final_shapes_supp.png', dpi=600, bbox_inches='tight')
    plt.savefig(out_dir + 'final_shapes_supp.pdf', bbox_inches='tight')
    print('wrote', out_dir + 'final_shapes_supp.png/.pdf')
