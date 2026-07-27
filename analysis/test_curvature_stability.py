#!/bin/python
"""
How stable is the extracted slab curvature K at 300 km?

Motivation.  The in-slab normal-stress term is K*Qn, and with the pressure
part included Qn = H*(sigma_n_bar - P_ext).  In the eta' = 1000 fixed-SP model
the slab is in extreme along-slab extension (tau_n_bar about -48 MPa) while
its curvature is small, so the term is a small, noisy number multiplied by a
very large stress.  That model is also the only one with a large residual in
the full force balance, confined to t >= 26 at 250 and 300 km.  If K is
poorly determined there, that explains the residual as a measurement limit
rather than as evidence against the pressure term.

Method.  K is extracted exactly as extract_properties.py does (llith = 0.5
contour -> midplane -> dip -> Savitzky-Golay smoothing -> centred difference
-> Savitzky-Golay again), but with the smoothing window swept.  The pipeline
hardcodes 601 in two places: the dip smoothing in extract_properties.py and
the K smoothing inside functions.get_curvature_slab_midplane.  Both are swept
together here, since both act on the same contour.

Because Qn does not depend on K, the full curvature term scales linearly with
K, so the implied change in the force-balance residual is
    delta_misfit = (K(window) - K(601)) * Qn
with Qn taken from the measured term at the reference window.

Usage:  python3 test_curvature_stability.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ANALYSIS = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, ANALYSIS)
from functions import create_grid, get_slab_midplane, get_dip_slab_midplane
from functions import get_curvature_at_certain_depth

c_llith_col = 25
x_col = 30
y_col = 31
xmax = 5800.e3
ymax = 1450.e3
analysis_depth = 300.e3
FIRST_TIME = 8

WINDOWS = [301, 401, 601, 801, 1001]
REF_WINDOW = 601

# (model, timesteps).  The suspect model plus two controls: a fixed-SP model
# with the same boundary condition but a weaker slab, and a free-plate model.
CASES = [
    ('2D_compositional_subd_lower-res_new_FixedSP_1000plates2', [24, 25, 26, 27, 29]),
    ('2D_compositional_subd_lower-res_new_FixedSP_250plates',   [24, 26, 29]),
    ('2D_compositional_subd_lower-res_new2',                    [24, 26, 29]),
]

# `python3 test_curvature_stability.py <model> <t> [<t> ...]` overrides the
# list above, so the sweep can be run over the whole suite in parallel.
if len(sys.argv) > 2:
    CASES = [(sys.argv[1], [int(a) for a in sys.argv[2:]])]


def curvature_with_window(llith_points, dips_unsmoothed, window):
    """Replicate the pipeline's K, with the smoothing window parameterised."""
    n = len(llith_points)
    w = min(window, n - 1 if (n - 1) % 2 else n - 2)   # savgol needs odd <= n
    if w < 5:
        return np.nan
    dips = savgol_filter(dips_unsmoothed[:, 0], w, 3)

    K_un = np.zeros((n, 1))
    for i in range(n):
        if i == 0:
            ds = llith_points[i + 1, 2] - llith_points[i, 2]
            ddip = dips[i + 1] - dips[i]
        elif i == n - 1:
            ds = llith_points[i, 2] - llith_points[i - 1, 2]
            ddip = dips[i] - dips[i - 1]
        else:
            ds = llith_points[i + 1, 2] - llith_points[i - 1, 2]
            ddip = dips[i + 1] - dips[i - 1]
        K_un[i] = np.deg2rad(ddip) / (ds * 1.e3)
    K = savgol_filter(K_un[:, 0], w, 3)
    K_mid, _, _, _ = get_curvature_at_certain_depth(K, K, llith_points,
                                                    analysis_depth)
    return K_mid


def main():
    xmin2, xmax2 = 2500.e3, 3750.e3
    X2, Y2 = create_grid(xmin2, xmax2, ymax - 600.e3, ymax, 0.25e3)

    print('CURVATURE STABILITY AT 300 km: K [1/km] vs Savitzky-Golay window')
    print('reference window = %d (the pipeline value)\n' % REF_WINDOW)

    for model, times in CASES:
        b = np.loadtxt(ANALYSIS + 'text_files/TESTB/%s.z300.0.shear-dz10.0.'
                       'ds10.0.prof-dz1.0km.txt' % model)
        tdat = np.atleast_2d(np.loadtxt(ANALYSIS + 'text_files/'
                             'curvature_pressure_term/%s.z300.0.txt' % model))
        tmap = {int(r[0]): r[15] for r in tdat}

        print('=' * 76)
        print(model.replace('2D_compositional_subd_', ''))
        print('%4s %5s' % ('t', 'npts') + ''.join('%11d' % w for w in WINDOWS)
              + '%12s %11s' % ('spread', 'd_misfit'))

        for t in times:
            csv = ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (model, t)
            if not os.path.exists(csv) or t not in tmap:
                print('%4d  (no data)' % t)
                continue
            md = pd.read_csv(csv).to_numpy()
            llith = griddata((md[:, x_col], md[:, y_col]), md[:, c_llith_col],
                             (X2, Y2), method='linear')
            ll = plt.contour(X2 / 1.e3, (ymax - Y2) / 1.e3, llith, levels=[0.5])
            pts = get_slab_midplane(ll.allsegs[0][0], 110., 575.)
            plt.close('all')
            dips_un = get_dip_slab_midplane(pts)

            Ks = [curvature_with_window(pts, dips_un, w) for w in WINDOWS]
            Kref = Ks[WINDOWS.index(REF_WINDOW)]

            row = b[t - FIRST_TIME]
            term_ref = (row[17] + tmap[t]) / 1.e6          # MPa, full K*Qn
            Qn = term_ref / Kref if Kref != 0 else np.nan  # MPa per (1/m)
            dmis = max(abs((k - Kref) * Qn) for k in Ks)

            print('%4d %5d' % (t, len(pts))
                  + ''.join('%11.5f' % (k * 1e3) for k in Ks)
                  + '%11.1f%% %10.2f' % (
                      100. * (max(Ks) - min(Ks)) / max(abs(Kref), 1e-30), dmis))
            del md
        print()


if __name__ == '__main__':
    main()
