#!/bin/python
"""
Does including the curvature-pressure term T improve the full force balance?

The top row of the DP-vs-DP compilation plots the total resisting force
    y = DP - dQs/ds + K*Qn        (columns 3, 6, 17 of the extraction files)
against slab buoyancy x = Bslab (column 4).  Equation 2 also contains
    T = (H/2)(dtheta/ds)(2*Pbar_n - P_subslab - P_wedge)
which the manuscript drops.  This script reports the misfit y - x with and
without T, at each analysis depth, so the effect on closure is quantified
rather than assumed.

Usage:  python3 compare_balance_with_curvature_pressure.py [depth_m ...]
        (default: 250e3 300e3 350e3)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TDIR = os.path.join(HERE, 'text_files', 'curvature_pressure_term')
BDIR = os.path.join(HERE, 'text_files', 'TESTB')

FIRST_TIME = 8
T_COL = 15          # T_cut
DP, BSLAB, DQDS, KQN = 3, 4, 6, 17

OVERTURNED = {
    '2D_compositional_subd_lower-res_new_FixedOP_375plates',
    '2D_compositional_subd_FixedOP_lower-res_new',
    '2D_compositional_subd_lower-res_new_1000plates',
    '2D_compositional_subd_lower-res_new_FixedOP_1000plates',
}

MODELS = [
    '2D_compositional_subd_lower-res_new_50plates',
    '2D_compositional_subd_lower-res_new_FixedSP_50plates',
    '2D_compositional_subd_lower-res_new_FixedOP_50plates',
    '2D_compositional_subd_lower-res_new_250plates',
    '2D_compositional_subd_lower-res_new_FixedSP_250plates',
    '2D_compositional_subd_lower-res_new_FixedOP_250plates',
    '2D_compositional_subd_lower-res_new_375plates',
    '2D_compositional_subd_lower-res_new_FixedSP_375plates',
    '2D_compositional_subd_lower-res_new_FixedOP_375plates',
    '2D_compositional_subd_lower-res_new2',
    '2D_compositional_subd_FixedSP_lower-res_new2',
    '2D_compositional_subd_FixedOP_lower-res_new',
    '2D_compositional_subd_lower-res_new_1000plates',
    '2D_compositional_subd_lower-res_new_FixedSP_1000plates2',
    '2D_compositional_subd_lower-res_new_FixedOP_1000plates',
]


def load_depth(depth):
    zkm = '%.1f' % (depth / 1.e3)
    old, new, ot, T = [], [], [], []
    for m in MODELS:
        tf = os.path.join(TDIR, '%s.z%s.txt' % (m, zkm))
        bf = os.path.join(BDIR, '%s.z%s.shear-dz10.0.ds10.0.prof-dz1.0km.txt'
                          % (m, zkm))
        if not os.path.exists(tf):
            sys.exit('missing %s (run many_curvature-pressure-term.sh %s)'
                     % (tf, depth))
        td = np.atleast_2d(np.loadtxt(tf))
        bd = np.loadtxt(bf)
        for r in td:
            b = bd[int(round(r[0])) - FIRST_TIME]
            o = (b[DP] - b[DQDS] + b[KQN] - b[BSLAB]) / 1e6
            t = r[T_COL] / 1e6
            old.append(o)
            new.append(o + t)
            T.append(t)
            ot.append(m in OVERTURNED)
    return (np.array(old), np.array(new), np.array(T), np.array(ot))


depths = [float(a) for a in sys.argv[1:]] or [250.e3, 300.e3, 350.e3]

print('=' * 78)
print('EFFECT OF THE CURVATURE-PRESSURE TERM ON FORCE-BALANCE CLOSURE')
print('misfit = (DP - dQs/ds + K*Qn [+ T]) - Bslab,  timesteps t >= 11')
print('=' * 78)

for depth in depths:
    old, new, T, ot = load_depth(depth)
    print('\n--- %.0f km ---' % (depth / 1.e3))
    for name, msk in (('all', np.ones(len(old), bool)),
                      ('non-overturned', ~ot)):
        o, n = old[msk], new[msk]
        rms_o, rms_n = np.sqrt((o ** 2).mean()), np.sqrt((n ** 2).mean())
        print('  %-15s N=%3d' % (name, msk.sum()))
        print('     as published : mean %+0.3f  median %+0.3f  RMS %.3f  '
              'max|.| %.3f MPa' % (o.mean(), np.median(o), rms_o, np.abs(o).max()))
        print('     including T  : mean %+0.3f  median %+0.3f  RMS %.3f  '
              'max|.| %.3f MPa' % (n.mean(), np.median(n), rms_n, np.abs(n).max()))
        print('     RMS change %+0.1f%%,  |misfit| reduced for %.0f%% of points,'
              '  median |T| = %.2f MPa'
              % (100. * (rms_n / rms_o - 1.), 100. * (np.abs(n) < np.abs(o)).mean(),
                 np.median(np.abs(T[msk]))))
print()
