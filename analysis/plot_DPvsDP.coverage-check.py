#!/bin/python
"""
Diagnostic: what the eta' = 1000 fixed-SP points look like under three options.

  (a) as published        column 17 = K*Qn with Qn = H*tau_n_bar
  (b) with T              column 17 = K*Qn with Qn = H*(sigma_n_bar - P_ext)
  (c) with T, covered     as (b), but timesteps where the slab midplane does
                          not resolve 300 km are dropped (column 38 of the
                          re-extraction, see many_property-extractions.recheck.sh)

300 km, the only depth the re-extraction covers.  Filled grey is the other
non-overturned models, hollow grey is the four overturned ones, red is
eta' = 1000 fixed-SP, and rings mark timesteps the coverage flag rejects.
Quoted statistics are for non-overturned models, as everywhere in the paper,
with the all-model values underneath.

Usage:  python3 plot_DPvsDP.coverage-check.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

ANALYSIS = os.path.dirname(os.path.abspath(__file__)) + '/'
SUF = '.z300.0.shear-dz10.0.ds10.0.prof-dz1.0km.txt'
TMIN = 3
SUSPECT = '2D_compositional_subd_lower-res_new_FixedSP_1000plates2'
OVERTURNED = ['FixedOP_375plates', 'FixedOP_lower-res_new',
              'lower-res_new_1000plates', 'FixedOP_1000plates']

font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
if os.path.exists(font_path):
    mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7


def force(row):
    """DP - dQs/ds + K*Qn, in MPa."""
    return (row[3] - row[6] + row[17]) / 1.e6


def main():
    import glob
    models = [os.path.basename(f)[:-len(SUF)]
              for f in sorted(glob.glob(ANALYSIS + 'text_files/TESTB' + '/*' + SUF))]
    isOT = lambda m: any(o in m for o in OVERTURNED)

    # coverage flag, per model per timestep, from the 300 km re-extraction
    flag = {}
    for m in models:
        p = ANALYSIS + 'text_files/TESTD_legacy/' + m + SUF
        if os.path.exists(p):
            d = np.loadtxt(p)
            flag[m] = {int(r[0]): r[38] for r in d}

    panels = [('a) as published', 'TESTB', False),
              ('b) with T', 'TESTB_withT', False),
              ('c) with T, covered only', 'TESTB_withT', True)]

    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.4))

    for ax, (title, src, cut) in zip(axes, panels):
        allv = []
        otv = []
        for m in models:
            d = np.loadtxt(ANALYSIS + 'text_files/' + src + '/' + m + SUF)[TMIN:]
            sus = (m == SUSPECT)
            ot = isOT(m)
            for r in d:
                t = int(r[0])
                ok = flag.get(m, {}).get(t, 1.0) == 1.0
                if cut and not ok:
                    continue
                x, y = r[4] / 1.e6, force(r)
                (otv if ot else allv).append(y - x)
                if ot:
                    ax.scatter(x, y, s=9, facecolors='none', edgecolor='0.62',
                               linewidth=0.5, zorder=3)
                else:
                    ax.scatter(x, y, s=14 if sus else 9,
                               color='firebrick' if sus else '0.72',
                               edgecolor='black' if (sus and not ok) else 'none',
                               linewidth=0.9 if (sus and not ok) else 0,
                               zorder=5 if sus else 2)

        ax.plot([-3, 33], [-3, 33], color='black', linewidth=1, zorder=1)
        ax.set_xlim(-3, 33)
        ax.set_ylim(-3, 33)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$B_{slab}$   [MPa]')
        ax.set_title(title, fontsize=8, loc='left')
        ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)

        v = np.array(allv)
        a = np.array(allv + otv)
        ax.annotate('non-overturned:  N = %d\nRMS = %.2f MPa,  max = %.2f MPa\n'
                    'all:  N = %d,  RMS = %.2f,  max = %.2f'
                    % (len(v), np.sqrt((v ** 2).mean()), np.abs(v).max(),
                       len(a), np.sqrt((a ** 2).mean()), np.abs(a).max()),
                    xy=(0.04, 0.96), xycoords='axes fraction', va='top',
                    ha='left', fontsize=6.5)

    axes[0].set_ylabel(r'$\Delta P  +  \sigma_{slab}$   [MPa]')
    axes[2].annotate('red: $\\eta$\' = 1000 fixed SP\n'
                     'hollow: overturned\nrings: 300 km not resolved',
                     xy=(0.96, 0.06), xycoords='axes fraction', va='bottom',
                     ha='right', fontsize=6.5)

    out = ANALYSIS + 'plots/DP-comparisons/compilations/DPvsDP.coverage-check.z300'
    plt.tight_layout()
    plt.savefig(out + '.png', bbox_inches='tight', dpi=400)
    plt.savefig(out + '.pdf', bbox_inches='tight')
    print('wrote ' + out + '.png')


if __name__ == '__main__':
    main()
