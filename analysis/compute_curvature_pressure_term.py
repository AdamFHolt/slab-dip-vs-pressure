#!/bin/python
"""
Compute the curvature-pressure term that Section 2.1 of the manuscript argues
away as negligible:

    T = (H/2) * (dtheta/ds) * (2*Pbar_n - P_subslab - P_wedge)

Pbar_n is the dynamic pressure averaged across the slab interior on the same
slab-normal cross-section that the existing pipeline uses to integrate the
in-slab deviatoric stresses.  The pipeline (extract_properties.py) already
interpolates this pressure profile as P_slab but never writes it out, so this
script reconstructs it from the archived csv_outputs.  No models are re-run.

Everything else is taken straight from the existing analysis:
  - dip, curvature K, slab-normal thickness H are read from the same
    text_files/TESTB/*.z300.0... files the paper's figures use,
  - P_subslab and P_wedge are re-extracted with get_nearslab_stresses using the
    same ds = 10 km stand-off, so all three pressures share one gauge (the
    combination 2*Pbar_n - P_subslab - P_wedge is gauge invariant, but the
    Pleft/Pright columns saved in the text files have a per-run reference
    pressure removed, so they cannot be mixed with a raw interior average),
  - the slab interior is the profile minus profcut = 10 km at each end, the
    same cut used for the deviatoric normal-stress resultant Q_n.

Usage:  python3 compute_curvature_pressure_term.py <model_name>
        (all 15 models in parallel via many_curvature-pressure-term.sh)
Writes: text_files/curvature_pressure_term/<model_name>.txt
"""
import os
import sys

import numpy as np
import pandas as pd
from numpy import trapz

ANALYSIS = '/home/holt/Projects/ASPECT/subd_2D/compositional/analysis/'
sys.path.insert(0, ANALYSIS)
from functions import extract_horiz_prof, get_slablocation_from_horiz_prof
from functions import get_nearslab_stresses, get_stress_profile

# csv column layout (as in extract_properties.py)
c_ulith_col = 24
c_llith_col = 25
P_col = 29          # "nonadiabatic_pressure" = dynamic pressure
x_col = 30
y_col = 31
sxx_col = 3
sxy_col = 4

# analysis settings, matching many_property-extractions.sh / all_plots.sh
ymax = 1450.e3
ds = 10.e3          # stand-off from slab surface for P_subslab / P_wedge
dz = 1.e3           # half-height of the horizontal profile slabs
profcut = 10.0      # km trimmed from each end of the cross-slab profile
first_time = 8      # first timestep in the extraction files
tactual_min = 11    # first timestep used by the paper's figures (tmin = 3)

TXT = 'text_files/TESTB/{model}.z{zkm}.shear-dz10.0.ds10.0.prof-dz1.0km.txt'
OUTDIR = ANALYSIS + 'text_files/curvature_pressure_term'

# text-file column indices
DIP_IND = 5
H_IND = 9
K_IND = 11
DP_IND = 3
KQN_IND = 17


def process_timestep(md, row, analysis_depth):
    """Return a dict of everything computed for one timestep."""
    dip = row[DIP_IND]
    H = row[H_IND]
    K = row[K_IND]

    prof = extract_horiz_prof(md, analysis_depth, ymax, dz, x_col, y_col)
    x_left, y_left, x_center, y_center, x_right, y_right = \
        get_slablocation_from_horiz_prof(prof, x_col, y_col, c_ulith_col,
                                         c_llith_col, ymax)

    sind = np.sin(np.deg2rad(dip))
    H_recon = (x_right - x_left) * sind
    thick_ll = (x_center - x_left) * sind
    thick_ul = (x_right - x_center) * sind

    # argument order follows extract_properties.py exactly
    nearslab = get_nearslab_stresses(
        y_center, thick_ul, thick_ll, dip, md, ds, ymax, x_col, y_col, P_col,
        c_llith_col, c_ulith_col, sxx_col, sxy_col, dz)
    Pleft, Pright = nearslab[0], nearslab[1]
    slab_x_left, slab_y_left, slab_x_right, slab_y_right = nearslab[10:]

    _, _, P_slab, profile = get_stress_profile(
        md, slab_x_left, slab_y_left, slab_x_right, slab_y_right, ymax,
        x_col, y_col, sxx_col, sxy_col, P_col)

    d = profile[:, 2]
    cut = np.where((d > profcut) & (d < d[-1] - profcut))[0]

    Pbar_cut = trapz(P_slab[cut], d[cut]) / (d[cut][-1] - d[cut][0])
    Pbar_full = trapz(P_slab, d) / (d[-1] - d[0])
    Pbar_med = np.median(P_slab[cut])

    bracket_cut = 2. * Pbar_cut - Pleft - Pright
    bracket_full = 2. * Pbar_full - Pleft - Pright
    bracket_med = 2. * Pbar_med - Pleft - Pright

    return dict(dip=dip, H=H, K=K, H_recon=H_recon,
                Pleft=Pleft, Pright=Pright, DP_recon=Pleft - Pright,
                Pbar_cut=Pbar_cut, Pbar_full=Pbar_full, Pbar_med=Pbar_med,
                bracket_cut=bracket_cut, bracket_full=bracket_full,
                bracket_med=bracket_med,
                T_cut=0.5 * H * K * bracket_cut,
                T_full=0.5 * H * K * bracket_full,
                T_med=0.5 * H * K * bracket_med,
                prof_span=d[-1])


def main():
    model = sys.argv[1]
    analysis_depth = float(sys.argv[2]) if len(sys.argv) > 2 else 300.e3
    zkm = '%.1f' % (analysis_depth / 1.e3)
    os.makedirs(OUTDIR, exist_ok=True)

    txt = np.loadtxt(ANALYSIS + TXT.format(model=model, zkm=zkm))
    n_rows = txt.shape[0]

    out = []
    for irow in range(tactual_min - first_time, n_rows):
        time = irow + first_time
        csv = ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (model, time)
        if not os.path.exists(csv):
            print('%s t=%d: MISSING CSV' % (model, time), flush=True)
            continue
        md = pd.read_csv(csv).to_numpy()
        r = process_timestep(md, txt[irow, :], analysis_depth)
        del md

        out.append([time, r['dip'], r['H'], r['K'], r['H_recon'],
                    r['Pleft'], r['Pright'], r['DP_recon'], txt[irow, DP_IND],
                    r['Pbar_cut'], r['Pbar_full'], r['Pbar_med'],
                    r['bracket_cut'], r['bracket_full'], r['bracket_med'],
                    r['T_cut'], r['T_full'], r['T_med'],
                    txt[irow, KQN_IND], r['prof_span']])
        print('%s t=%2d  HK/2=%7.4f  bracket=%8.3f MPa  T=%8.4f MPa  '
              '(dDP=%.2e MPa, dH=%.2e km)'
              % (model, time, 0.5 * r['H'] * r['K'], r['bracket_cut'] / 1e6,
                 r['T_cut'] / 1e6,
                 abs(r['DP_recon'] - txt[irow, DP_IND]) / 1e6,
                 abs(r['H_recon'] - r['H']) / 1e3), flush=True)

    header = ('time dip H K H_recon Pleft Pright DP_recon DP_txt '
              'Pbar_cut Pbar_full Pbar_med bracket_cut bracket_full '
              'bracket_med T_cut T_full T_med KQn prof_span_km')
    np.savetxt(os.path.join(OUTDIR, '%s.z%s.txt' % (model, zkm)),
               np.array(out), header=header)


if __name__ == '__main__':
    main()
