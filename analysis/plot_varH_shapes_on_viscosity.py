#!/usr/bin/env python3
"""
Fig.-2 grid (plate configuration x eta' = 50/250/375/500/1000) of the
H = 100 km viscosity snapshots, with the H = 70 and H = 130 km slab outlines
drawn on top of the three eta' that have them.

Each panel is the H = 100 model at the output Fig. 2 uses and the same
1200 x 1000 km zoom (many_simple-plots_viscosity.sh), lajolla_r over
log10(eta) = 21.5-24, Fig. 2 fonts.  The outlines are the 0.5 contour of the
whole lithosphere composition (crust + upper + lower lithosphere), taken at
the output whose slab-tip depth is closest to the H = 100 tip in the snapshot
(text_files/slab_tip_depths.txt) and shifted horizontally so the slab top at
250 km depth coincides with the H = 100 slab.

Usage: python3 plot_varH_shapes_on_viscosity.py   (from analysis/)
Writes plots/evolution/simple/varH_shapes_on_viscosity.png/.pdf (3x5 composite) and one
PDF/PNG per panel in plots/evolution/simple/zoomed/varH_overlay/.
"""
import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from cmcrameri import cm as cmc
from scipy.interpolate import griddata
from multiprocessing import Pool

# ---- style, as plot_onestep_simple-viscosity.zoomed.py ----
font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
if os.path.exists(font_path):
    fm.fontManager.addfont(font_path)
    mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams.update({'font.size': 7, 'axes.labelsize': 7, 'axes.labelpad': 1.25,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'xtick.major.pad': 2, 'ytick.major.pad': 2,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.25, 'ytick.minor.size': 1.25})

P = '2D_compositional_subd_lower-res_new_'
YMAX = 1450e3
ZREF = 250.          # km, alignment depth
ETAS = [50, 250, 375, 500, 1000]; BCS = ['free', 'FixedSP', 'FixedOP']
# H = 100 panels exactly as in many_simple-plots_viscosity.sh: (model, output, x_center_km)
H100 = {('free',50): (P+'50plates', 16, 3250), ('free',250): (P+'250plates', 16, 3250), ('free',375): (P+'375plates', 17, 3250),
        ('free',500): ('2D_compositional_subd_lower-res_new2', 17, 3350), ('free',1000): (P+'1000plates', 18, 3250),
        ('FixedSP',50): (P+'FixedSP_50plates', 16, 3200), ('FixedSP',250): (P+'FixedSP_250plates', 16, 3000), ('FixedSP',375): (P+'FixedSP_375plates', 17, 3000),
        ('FixedSP',500): ('2D_compositional_subd_FixedSP_lower-res_new2', 17, 2900), ('FixedSP',1000): (P+'FixedSP_1000plates2', 18, 2900),
        ('FixedOP',50): (P+'FixedOP_50plates', 16, 3250), ('FixedOP',250): (P+'FixedOP_250plates', 16, 3350), ('FixedOP',375): (P+'FixedOP_375plates', 17, 3350),
        ('FixedOP',500): ('2D_compositional_subd_FixedOP_lower-res_new', 17, 3350), ('FixedOP',1000): (P+'FixedOP_1000plates', 18, 3350)}
def varH(bc, eta, H): return P + ('' if bc == 'free' else bc + '_') + '%dplates_H%dkm' % (eta, H)
TICK_FS, LABEL_FS = 6, 7          # Fig. 2 sizes, used for the per-panel files
STYLE = {70: dict(color='deepskyblue', lw=0.9), 130: dict(color='limegreen', lw=0.9)}

tips = {}
for line in open('text_files/slab_tip_depths.txt'):
    p = line.split()
    if p and not p[0].startswith('#'): tips[(p[0], int(p[1]))] = float(p[2])

import re
def model_time_myr(m, t):
    pvd = 'raw_outputs/%s/solution.pvd' % m
    if os.path.exists(pvd):
        return [float(x) for x in re.findall(r'timestep="([^"]+)"', open(pvd).read())][t] / 1e6
    # as plot_onestep_simple-viscosity.zoomed.py: statistics line 16 + 50 per output
    return float(open('raw_outputs/%s/statistics' % m).readlines()[16 + 50 * t].split()[1]) / 1e6

def matched_output(m, target):
    cand = [(abs(t - target), i) for (mm, i), t in tips.items() if mm == m]
    return min(cand)[1] if cand else None

def load(m, t):
    fn = 'csv_outputs/%s/full.%d.csv' % (m, t)
    hdr = open(fn).readline().strip().replace('"', '').split(',')
    d = np.loadtxt(fn, delimiter=',', skiprows=1)
    cols = [hdr.index(c) for c in ('Points:0', 'Points:1', 'crust', 'ulith', 'llith', 'viscosity')]
    return d, cols

def outline(d, cols, x0, x1):
    xc, yc, cc, uc, lc, _ = cols
    lith = d[:, cc] + d[:, uc] + d[:, lc]
    X, Y = np.meshgrid(np.arange(x0 - 400e3, x1 + 400e3, 2e3), np.arange(YMAX - 1100e3, YMAX + 1, 2e3))
    ll = griddata((d[:, xc], d[:, yc]), lith, (X, Y), method='linear')
    cs = plt.contour(X / 1e3, (YMAX - Y) / 1e3, ll, levels=[0.5]); plt.close('all')
    seg = max(cs.allsegs[0], key=len)
    near = seg[np.abs(seg[:, 1] - ZREF) < 4]
    return seg, near[:, 0].max()          # slab top (OP side) at ZREF

def panel(key):
    bc, eta = key
    m100, t, xcen = H100[key]
    x0, x1 = (xcen - 600) * 1e3, (xcen + 600) * 1e3
    d, cols = load(m100, t)
    X, Y = np.meshgrid(np.arange(x0, x1 + 1, 2e3), np.arange(YMAX - 1000e3, YMAX + 1, 2e3))
    visc = griddata((d[:, cols[0]], d[:, cols[1]]), d[:, cols[5]], (X, Y), method='nearest')
    out = {}
    if eta in (250, 500, 1000):
        _, xref100 = outline(d, cols, x0, x1)
        target = tips[(m100, t)]
        for H in (70, 130):
            m = varH(bc, eta, H); tm = matched_output(m, target)
            if tm is None or not os.path.exists('csv_outputs/%s/full.%d.csv' % (m, tm)): continue
            dd, cc = load(m, tm)
            seg, xref = outline(dd, cc, x0, x1)
            out[H] = (seg, xref100 - xref, tm, tips[(m, tm)])
    return key, X / 1e3, (YMAX - Y) / 1e3, np.log10(visc), out, (x0 / 1e3, x1 / 1e3), t, tips.get((m100, t), np.nan), model_time_myr(m100, t)

def draw(ax, X, Z, lv, out, xlim):
    ax.contourf(X, Z, lv, cmap=cmc.lajolla_r, levels=np.linspace(21.5, 24, 501), extend='both', zorder=-1); ax.set_rasterization_zorder(0)   # viscosity field as an image in the PDFs, outlines stay vector
    for H in (70, 130):
        if H in out:
            seg, shift, _, _ = out[H]
            ax.plot(seg[:, 0] + shift, seg[:, 1], **STYLE[H])
    ax.set_xlim(xlim); ax.set_ylim(1000, 0); ax.set_aspect('equal', adjustable='box')
    ax.tick_params(direction='out', length=2, labelsize=TICK_FS)
    ax.set_xlabel('x [km]', fontsize=LABEL_FS); ax.set_ylabel('depth [km]', fontsize=LABEL_FS)

if __name__ == '__main__':
    with Pool(15) as pool:
        res = pool.map(panel, list(H100))
    os.makedirs('plots/evolution/simple/zoomed/varH_overlay', exist_ok=True)
    TICK_FS, LABEL_FS = 7, 8          # a touch larger in the composite
    fig, axs = plt.subplots(3, 5, figsize=(11, 5.4))
    BC_LABEL = {'free': 'free plates', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
    for key, X, Z, lv, out, xlim, t, tip, myr in res:
        bc, eta = key; i, j = BCS.index(bc), ETAS.index(eta); ax = axs[i, j]
        draw(ax, X, Z, lv, out, xlim)
        ax.annotate('%.1f Myr' % myr, xy=(0.03, 0.04), xycoords='axes fraction', ha='left', va='bottom', fontsize=8)
        if i == 0:
            ax.set_title("$\\eta'$ = %d" % eta, fontsize=11, pad=6)
        if j == 0:
            ax.annotate(BC_LABEL[bc], xy=(-0.42, 0.5), xycoords='axes fraction', ha='center', va='center',
                        rotation=90, fontsize=11)
        if i < 2:
            ax.set_xlabel(''); ax.tick_params(labelbottom=False)
        if j > 0:
            ax.set_ylabel(''); ax.tick_params(labelleft=False)
        if i == 0 and j == 1:
            from matplotlib.lines import Line2D
            leg = ax.legend(handles=[Line2D([], [], label='70 km', **STYLE[70]), Line2D([], [], label='130 km', **STYLE[130])],
                            title='Variable H', loc='lower left', bbox_to_anchor=(0.03, 0.36), fontsize=7, title_fontsize=7.5, frameon=True, facecolor='none',
                            edgecolor='black', fancybox=False, handlelength=1.4, borderpad=0.3, labelspacing=0.2, handletextpad=0.5, borderaxespad=0.25)
            leg.get_frame().set_linewidth(0.6); leg.get_frame().set_edgecolor('black')
        line = '%-8s %5d: H100 out %d tip %.0f' % (bc, eta, t, tip)
        for H in (70, 130):
            if H in out: line += ' | H%d out %d tip %.0f shift %+.0f km' % (H, out[H][2], out[H][3], out[H][1])
        print(line + ' | %.1f Myr' % myr)
        TICK_FS, LABEL_FS = 6, 7
        f1, a1 = plt.subplots(figsize=(3.4, 3.4))
        draw(a1, X, Z, lv, out, xlim)
        TICK_FS, LABEL_FS = 7, 8
        base = 'plots/evolution/simple/zoomed/varH_overlay/%s.%d.varH' % (H100[key][0], H100[key][1])
        f1.savefig(base + '.png', bbox_inches='tight', dpi=500); f1.savefig(base + '.pdf', bbox_inches='tight', dpi=400); plt.close(f1)
    fig.tight_layout(w_pad=0.3, h_pad=0.0); fig.savefig('plots/evolution/simple/varH_shapes_on_viscosity.png', dpi=200, bbox_inches='tight'); fig.savefig('plots/evolution/simple/varH_shapes_on_viscosity.pdf', bbox_inches='tight', dpi=400)
    print('wrote plots/evolution/simple/varH_shapes_on_viscosity.png/.pdf and per-panel files in plots/evolution/simple/zoomed/varH_overlay/')
