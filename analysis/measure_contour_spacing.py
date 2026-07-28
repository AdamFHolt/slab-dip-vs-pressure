#!/bin/python
"""
How much ARC LENGTH does the 601-point Savitzky-Golay window actually span?

The pipeline smooths dip and curvature with a window of a fixed number of
contour points (601, hardcoded in extract_properties.py and in
functions.get_curvature_slab_midplane).  The midplane contour is resampled by
matplotlib's contour generator, so its point spacing is not constant across
models or through time: when a slab necks, the same depth range is described by
fewer points, each one further apart, and the window quietly widens in physical
terms.  That is the mechanism behind the eta' = 1000 fixed-SP curvature
instability (see test_curvature_stability.py).

This script measures, for every model and a few timesteps, the number of
midplane points, the arc length they span, the median point spacing, and hence
the physical length of a 601-point window.  The spread in that last number is
what an arc-length window is meant to remove.

Usage:  python3 measure_contour_spacing.py <model> [<t> ...]
        (writes one line per timestep to stdout)
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ANALYSIS = os.path.dirname(os.path.abspath(__file__)) + '/'
sys.path.insert(0, ANALYSIS)
from functions import create_grid, get_slab_midplane

c_llith_col = 25
x_col = 30
y_col = 31
ymax = 1450.e3
REF_WINDOW = 601


def main():
    model = sys.argv[1]
    times = [int(a) for a in sys.argv[2:]]

    X2, Y2 = create_grid(2500.e3, 3750.e3, ymax - 600.e3, ymax, 0.25e3)

    for t in times:
        csv = ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (model, t)
        if not os.path.exists(csv):
            continue
        md = pd.read_csv(csv).to_numpy()
        llith = griddata((md[:, x_col], md[:, y_col]), md[:, c_llith_col],
                         (X2, Y2), method='linear')
        cont = plt.contour(X2 / 1.e3, (ymax - Y2) / 1.e3, llith, levels=[0.5])
        pts = get_slab_midplane(cont.allsegs[0][0], 110., 575.)
        plt.close('all')
        del md

        s = pts[:, 2]                       # cumulative arc length, km
        n = len(s)
        span = s[-1] - s[0]
        dsm = np.median(np.diff(s))
        # matplotlib's own spacing, i.e. what 601 points really covers
        print('%-55s %3d %6d %8.1f %8.4f %8.1f'
              % (model.replace('2D_compositional_subd_', ''), t, n, span,
                 dsm, REF_WINDOW * dsm))
        sys.stdout.flush()


if __name__ == '__main__':
    main()
