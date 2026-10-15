#!/bin/python
"""
Fig. 4 for the revision: the time-averaged force balance of the H = 100 km
models, split into its parts, against eta' for each plate configuration.
Replaces plot_forces-bars.no-norm.py (kept to reproduce the submitted figure).

Each bar stacks the time-averaged terms, positive ones up from zero and
negative ones down: dQ/ds (dark blue) and K*Qn (light blue), the two parts of
sigma_slab, then DP (red).  Black markers are the averaged model total
DP + sigma_slab with its standard deviation over time (in MPa; the submitted
figure plotted the stdev of total/B_slab on the MPa axis); grey markers are
B_slab.  Marker shape is the plate configuration as in Fig. 3.  Overturned
slabs (the set in plot_forces.varH.py, i.e. slab-top dip > 90 deg at 660 km)
have hatched bars and hollow markers.  The strip on top is the mean signed K.

Data and averaging window are those of Fig. 3 (plot_forces.varH.py closure):
text_files/withT_varH/ (post window-fix extraction), outputs from 11 to the
end of the extraction window, coverage-flagged rows dropped, the fixed-SP
eta' = 50 run trimmed to output 22.  Prints the averages for the caption.

Usage:  python3 plot_forces-bars.fig4.py 300.e3 10.e3 10.e3 1.e3
Writes  plots/DP-comparisons/compilations/manuscript_figs/forces-bars.H100<suffix>.tmin3.png/pdf
"""
import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams.update({'font.size': 7, 'axes.labelsize': 10, 'xtick.labelsize': 8.5,
                     'ytick.labelsize': 8.5, 'hatch.linewidth': 0.6})

analysis_depth    = float(sys.argv[1])
analysis_depth_dz = float(sys.argv[2])
ds                = float(sys.argv[3])
dz                = float(sys.argv[4])

tactual_min = 11
tmin = tactual_min - 8
suffix = '.z%s.shear-dz%s.ds%s.prof-dz%skm' % (analysis_depth/1.e3, analysis_depth_dz/1.e3, ds/1.e3, dz/1.e3)
out_dir = 'plots/DP-comparisons/compilations/manuscript_figs/'
os.makedirs(out_dir, exist_ok=True)
# "varH": the supplement version, the H = 70 and 130 km models in two stacked
# panels like b) (eta' = 250, 500, 1000 only; no curvature strip)
varH_fig = 'varH' in sys.argv[5:]
plot_name = out_dir + ('forces-bars.H70-130' if varH_fig else 'forces-bars.H100') + suffix + '.tmin%d' % tmin
YLIM_VARH, YTICKS_VARH = (-5, 40), [0, 10, 20, 30, 40]

K_ind = 11; ss_ind = 17; sn_ind = 6; anal_ind = 4; DP_ind = 3; cov_ind = 38

# as plot_forces.varH.py
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
# the OVERTURNED set of plot_forces.varH.py (Fig. 3), keyed (plates, eta', H)
OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130),
              ('FixedOP', 375, 100), ('FixedOP', 500, 100), ('FixedOP', 1000, 100),
              ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70)}
LAST_OUTPUT = {('FixedSP', 50, 100): 22}

BCS  = ['FixedSP', 'free', 'FixedOP']          # left to right within each eta' group
ETAS = [250, 500, 1000] if varH_fig else [50, 250, 375, 500, 1000]
HS   = [70, 130] if varH_fig else [100]
bc_label  = {'free': 'free plates', 'FixedSP': 'fixed SP', 'FixedOP': 'fixed OP'}
bc_marker = {'FixedSP': 'v', 'free': 'o', 'FixedOP': '^'}
# x is eta' itself, except that 1000 is drawn at 725 (as in the submitted
# figure) with an axis break between 500 and 1000
XPOS   = {50: 50, 250: 250, 375: 375, 500: 500, 1000: 725}
XBREAK = 625
XLIM   = (-25, 810)
if varH_fig:
    # three evenly spaced clusters (no true scale, so no break); the figure is
    # just wide enough for the one-row legend and the overturned key beside it
    XPOS, XBREAK, XLIM = {250: 0, 500: 150, 1000: 300}, None, (-75, 375)
offset = {'FixedSP': -34., 'free': 0., 'FixedOP': 34.}
w = 24.
c_DP, c_dqds, c_KQn = 'red', 'mediumblue', 'cornflowerblue'

def averages(bc, eta, H):
    d = np.loadtxt('text_files/withT_varH/' + name(bc, eta, H) + suffix + '.txt')[tmin:]
    d = d[d[:, cov_ind] == 1]
    d = d[d[:, 0] <= LAST_OUTPUT.get((bc, eta, H), 1e9)]
    DP    = d[:, DP_ind] / 1.e6
    dqds  = -d[:, sn_ind] / 1.e6               # paper sign
    KQn   = d[:, ss_ind] / 1.e6
    total = DP + dqds + KQn
    return dict(N=len(d), DP=DP.mean(), dqds=dqds.mean(), KQn=KQn.mean(), total=total.mean(),
                std=total.std(), B=d[:, anal_ind].mean() / 1.e6, K=d[:, K_ind].mean() * 1.e6)

avg = {(bc, e, H): averages(bc, e, H) for H in HS for bc in BCS for e in ETAS}

print('%-8s %5s %4s %3s %7s %7s %7s %7s %6s %7s %8s' % ('plates', 'eta', 'H', 'N', 'DP', 'dQ/ds', 'KQn', 'total', 'std', 'B_slab', 'K e-3/km'))
for (bc, e, H), a in avg.items():
    print('%-8s %5d %4d %3d %7.2f %7.2f %7.2f %7.2f %6.2f %7.2f %8.3f%s'
          % (bc, e, H, a['N'], a['DP'], a['dqds'], a['KQn'], a['total'], a['std'], a['B'], a['K'],
             '  overturned' if (bc, e, H) in OVERTURNED else ''))

def bars(ax, axK, H):
    for e in ETAS:
        for bc in BCS:
            a = avg[(bc, e, H)]; x = XPOS[e] + offset[bc]; ot = (bc, e, H) in OVERTURNED     # hatched bars; markers as for the rest
            hatch = dict(hatch='/////', edgecolor='white', linewidth=0) if ot else dict(linewidth=0)
            up, down = 0., 0.
            for key, c in (('dqds', c_dqds), ('KQn', c_KQn), ('DP', c_DP)):
                v = a[key]
                if v >= 0:
                    ax.bar(x, v, bottom=up, width=w, color=c, zorder=3, **hatch); up += v
                else:
                    ax.bar(x, v, bottom=down, width=w, color=c, zorder=3, **hatch); down += v
            m = bc_marker[bc]
            ax.plot([x, x], [a['total'] - a['std'], a['total'] + a['std']], color='silver', linewidth=1, zorder=5)
            ax.scatter(x, a['B'], s=42, color='grey', edgecolor='grey', linewidth=0.25, marker=m, zorder=5)
            ax.scatter(x, a['total'], s=42, color='black', edgecolor='black', linewidth=0.25, marker=m, zorder=6)
            if axK is not None:
                axK.scatter(x, a['K'], s=36, color='black', edgecolor='black', linewidth=0.25, marker=m, zorder=6)

def frame(a_, ref):
    a_.axhline(0, color='silver', linestyle='--', linewidth=1, zorder=1)
    a_.grid(True, axis='y', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
    a_.set_xlim(*XLIM)
    if XBREAK is None:
        return
    # axis break: a white gap in the bottom and top spines with two short slashes,
    # the same size in every panel (height in axes units scales with 1/axes height)
    hy = 0.035 * ref.get_position().height / a_.get_position().height
    for yb in (0, 1):
        tr = mpl.transforms.blended_transform_factory(a_.transData, a_.transAxes)
        a_.add_patch(mpl.patches.Rectangle((XBREAK - 6, yb - hy), 12, 2 * hy, transform=tr,
                                           facecolor='white', edgecolor='none', clip_on=False, zorder=20))
        for xs in (XBREAK - 6, XBREAK + 6):
            a_.plot([xs - 5, xs + 5], [yb - hy, yb + hy], transform=tr, color='black',
                    linewidth=0.8, clip_on=False, zorder=21)

def stress_axis(ax, ylim, yticks):
    ax.set_xticks([XPOS[e] for e in ETAS]); ax.set_xticklabels([str(e) for e in ETAS])
    ax.set_ylabel(r'$\sigma$   [MPa]')
    ax.set_ylim(*ylim); ax.set_yticks(yticks)
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))

def legends(ax, ncol):
    # one boxed legend in one row across the top (forces only), clear of the bars;
    # the marker shapes (plate configuration) are explained in the caption
    terms = [Patch(facecolor=c_DP, label=r'$\Delta P$'),
             Patch(facecolor=c_dqds, label=r'$dQ/ds$'),
             Patch(facecolor=c_KQn, label=r'$K Q_n$')]
    marks = [Line2D([], [], marker='o', color='black', linestyle='', markersize=5.5, label=r'$\Delta P + \sigma_{slab}$'),
             Line2D([], [], marker='o', color='grey', linestyle='', markersize=5.5, label=r'$B_{slab}$')]
    leg = ax.legend(handles=terms + marks, ncol=ncol, loc='upper left', fontsize=8, frameon=True,
                    fancybox=False, edgecolor='black', facecolor='white', framealpha=1, borderpad=0.5,
                    labelspacing=0.35, handletextpad=0.4, handlelength=1.2, columnspacing=0.8,
                    borderaxespad=0.6)
    leg.get_frame().set_linewidth(0.5)
    leg.set_zorder(10)
    # the overturned key on its own, just outside the box on its right, centred on it
    ax.add_artist(leg)
    ax.figure.canvas.draw()
    bb = leg.get_window_extent().transformed(ax.transAxes.inverted())
    ax.legend(handles=[Patch(facecolor='gray', hatch='/////', edgecolor='white', linewidth=0, label='overturned')],
              loc='center left', bbox_to_anchor=(bb.x1 + 0.005, 0.5 * (bb.y0 + bb.y1)), bbox_transform=ax.transAxes,
              fontsize=8, frameon=False, borderaxespad=0, handletextpad=0.4, handlelength=1.2, handleheight=1.2)

if varH_fig:
    # supplement: the H = 70 and 130 km models, one panel each, no curvature strip
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 5.6), sharex=True, gridspec_kw=dict(hspace=0.2))
    for ax, H, lab in zip(axs, HS, 'ab'):
        bars(ax, None, H)
        frame(ax, axs[0])
        stress_axis(ax, YLIM_VARH, YTICKS_VARH)
        ax.set_title(lab + ')', loc='left', fontsize=10)
        ax.set_title('H = %d km' % H, fontsize=10)
    axs[1].set_xlabel(r"$\eta'$  ($\eta_{slab}$ / $\eta_{mantle}$)")
    legends(axs[0], 5)
else:
    fig, (axK, ax) = plt.subplots(2, 1, figsize=(7.0, 3.6), sharex=True,
                                  gridspec_kw=dict(height_ratios=[1, 4], hspace=0.22))
    bars(ax, axK, 100)
    for a_ in (ax, axK):
        frame(a_, ax)
    stress_axis(ax, (-5, 28), [-5, 0, 5, 10, 15, 20, 25])
    ax.set_xlabel(r"$\eta'$  ($\eta_{slab}$ / $\eta_{mantle}$)")
    axK.set_ylabel('$K$\n[10$^{-3}$ km$^{-1}$]')
    axK.set_ylim(-0.5, 3.0); axK.set_yticks([0, 1, 2, 3])
    axK.tick_params(axis='x', length=0)
    axK.set_title('a)', loc='left', fontsize=10); ax.set_title('b)', loc='left', fontsize=10)
    legends(ax, 5)

plt.savefig(plot_name + '.png', dpi=300, bbox_inches='tight')
plt.savefig(plot_name + '.pdf', bbox_inches='tight')
print('wrote', plot_name + '.png/.pdf')
