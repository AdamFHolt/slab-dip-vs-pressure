#!/bin/python
"""
Summarize the curvature-pressure term T computed by
compute_curvature_pressure_term.py.

T = (H/2) * (dtheta/ds) * (2*Pbar_n - P_subslab - P_wedge)

Reported against the slab buoyancy per unit slab area, B = drho*g*H, which is
49.05 MPa in the models (drho = 50 kg/m3, g = 9.81, H = 100 km), the same
normalization used for Lambda in the manuscript.

Usage:  python3 summarize_curvature_pressure_term.py [analysis_depth_m]
        (default 300e3, the depth used throughout the manuscript)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
INDIR = os.path.join(HERE, 'text_files', 'curvature_pressure_term')

ANALYSIS_DEPTH = float(sys.argv[1]) if len(sys.argv) > 1 else 300.e3
ZKM = '%.1f' % (ANALYSIS_DEPTH / 1.e3)
ZSUF = '.z%s.txt' % ZKM

B_NORM = 50. * 9.81 * 100.e3 * 1e-6   # 49.05 MPa

MODELS = [
    ('2D_compositional_subd_lower-res_new_50plates', 50, 'bothfree', False),
    ('2D_compositional_subd_lower-res_new_FixedSP_50plates', 50, 'fixedSP', False),
    ('2D_compositional_subd_lower-res_new_FixedOP_50plates', 50, 'fixedOP', False),
    ('2D_compositional_subd_lower-res_new_250plates', 250, 'bothfree', False),
    ('2D_compositional_subd_lower-res_new_FixedSP_250plates', 250, 'fixedSP', False),
    ('2D_compositional_subd_lower-res_new_FixedOP_250plates', 250, 'fixedOP', False),
    ('2D_compositional_subd_lower-res_new_375plates', 375, 'bothfree', False),
    ('2D_compositional_subd_lower-res_new_FixedSP_375plates', 375, 'fixedSP', False),
    ('2D_compositional_subd_lower-res_new_FixedOP_375plates', 375, 'fixedOP', True),
    ('2D_compositional_subd_lower-res_new2', 500, 'bothfree', False),
    ('2D_compositional_subd_FixedSP_lower-res_new2', 500, 'fixedSP', False),
    ('2D_compositional_subd_FixedOP_lower-res_new', 500, 'fixedOP', True),
    ('2D_compositional_subd_lower-res_new_1000plates', 1000, 'bothfree', True),
    ('2D_compositional_subd_lower-res_new_FixedSP_1000plates2', 1000, 'fixedSP', False),
    ('2D_compositional_subd_lower-res_new_FixedOP_1000plates', 1000, 'fixedOP', True),
]

# column indices written by compute_curvature_pressure_term.py
C = dict(time=0, dip=1, H=2, K=3, H_recon=4, Pleft=5, Pright=6, DP_recon=7,
         DP_txt=8, Pbar_cut=9, Pbar_full=10, Pbar_med=11, bracket_cut=12,
         bracket_full=13, bracket_med=14, T_cut=15, T_full=16, T_med=17,
         KQn=18, span=19)


def stats(x, label, unit=''):
    x = np.asarray(x)
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    print('  %-28s N=%3d  median=%8.4f  IQR=[%.4f, %.4f]  max=%8.4f %s'
          % (label, len(x), med, q1, q3, x.max(), unit))
    return med, q1, q3, x.max()


rows, visc, bc, ot = [], [], [], []
for model, eta, cond, overturned in MODELS:
    f = os.path.join(INDIR, model + ZSUF)
    if not os.path.exists(f):
        print('MISSING: %s' % model)
        continue
    d = np.atleast_2d(np.loadtxt(f))
    rows.append(d)
    visc.append(np.full(d.shape[0], eta))
    bc.append(np.full(d.shape[0], cond))
    ot.append(np.full(d.shape[0], overturned))

d = np.vstack(rows)
visc = np.concatenate(visc)
bc = np.concatenate(bc)
ot = np.concatenate(ot)

T = d[:, C['T_cut']] / 1e6            # MPa
T_full = d[:, C['T_full']] / 1e6
T_med = d[:, C['T_med']] / 1e6
bracket = d[:, C['bracket_cut']] / 1e6
HK2 = 0.5 * d[:, C['H']] * d[:, C['K']]
KQn = d[:, C['KQn']] / 1e6
B_actual = 50. * 9.81 * d[:, C['H']] * 1e-6

print('=' * 78)
print('CURVATURE-PRESSURE TERM  T = (H/2)(dtheta/ds)(2*Pbar_n - P_sub - P_wedge)')
print('%s km depth, timesteps t >= 11, B = %.2f MPa' % (ZKM, B_NORM))
print('=' * 78)

print('\nSANITY CHECKS (reconstruction vs archived extraction)')
print('  max |H_recon - H|          = %.3e km'
      % (np.abs(d[:, C['H_recon']] - d[:, C['H']]).max() / 1e3))
print('  max |DP_recon - DP_txt|    = %.3e MPa'
      % (np.abs(d[:, C['DP_recon']] - d[:, C['DP_txt']]).max() / 1e6))
print('  cross-slab profile span    = %.1f to %.1f km'
      % (d[:, C['span']].min(), d[:, C['span']].max()))

for name, mask in (('ALL TIMESTEPS', np.ones(len(T), bool)),
                   ('NON-OVERTURNED ONLY', ~ot),
                   ('OVERTURNED ONLY', ot)):
    if mask.sum() == 0:
        continue
    print('\n%s  (N = %d, %d models)' % (name, mask.sum(), len(set(visc[mask].astype(str) + bc[mask]))))
    stats(np.abs(T[mask]) / B_NORM, '|T| / B', '')
    stats(np.abs(T[mask]), '|T|', 'MPa')
    stats(np.abs(HK2[mask]), '(H/2)|dtheta/ds|', '')
    stats(np.abs(bracket[mask]), '|bracket|', 'MPa')
    stats(np.abs(KQn[mask]), '|K*Qn| (existing term)', 'MPa')
    r = np.abs(T[mask]) / B_NORM
    print('    fraction |T|/B > 0.01 : %5.1f%%' % (100. * (r > 0.01).mean()))
    print('    fraction |T|/B > 0.05 : %5.1f%%' % (100. * (r > 0.05).mean()))
    print('    fraction |T|/B > 0.10 : %5.1f%%' % (100. * (r > 0.10).mean()))
    print('    signed T: mean %+0.4f MPa, range %+0.4f to %+0.4f MPa'
          % (T[mask].mean(), T[mask].min(), T[mask].max()))

print('\nRELATION TO THE RETAINED DEVIATORIC TERM K*Qn')
print('  Traction continuity across the slab faces plus 2-D incompressibility')
print('  predicts Pbar_n - P_avg = tau_n, i.e. T should equal K*Qn (making the')
print('  full curvature contribution 2*K*Qn rather than K*Qn).')
for name, mask in (('all timesteps', np.ones(len(T), bool)),
                   ('non-overturned', ~ot)):
    good = mask & (np.abs(KQn) > 1e-3)
    ratio = T[good] / KQn[good]
    print('  %-16s N=%3d  median T/(K*Qn) = %+0.3f  IQR=[%+0.3f, %+0.3f]  '
          'corr(T, K*Qn) = %+0.3f'
          % (name, good.sum(), np.median(ratio),
             *np.percentile(ratio, [25, 75]),
             np.corrcoef(T[good], KQn[good])[0, 1]))

print('\nROBUSTNESS OF THE INTERIOR AVERAGE (non-overturned)')
m = ~ot
for lab, arr in (('interior cut (10 km trim)', T), ('full profile', T_full),
                 ('median instead of mean', T_med)):
    print('  %-28s median |T|/B = %.4f   max |T|/B = %.4f'
          % (lab, np.median(np.abs(arr[m])) / B_NORM,
             np.abs(arr[m]).max() / B_NORM))
print('  using per-timestep B = drho*g*H_actual (median H = %.1f km):'
      % (np.median(d[m, C['H']]) / 1e3))
print('      median |T|/B = %.4f   max |T|/B = %.4f'
      % (np.median(np.abs(T[m]) / B_actual[m]), (np.abs(T[m]) / B_actual[m]).max()))

print('\nPER-MODEL BREAKDOWN (median values)')
print('  %-10s %-9s %-4s %4s %10s %10s %10s %10s'
      % ('eta_prime', 'condition', 'ot', 'N', 'med|T|/B', 'med|T|MPa',
         'med HK/2', 'med|K*Qn|'))
for model, eta, cond, overturned in MODELS:
    f = os.path.join(INDIR, model + ZSUF)
    if not os.path.exists(f):
        continue
    dd = np.atleast_2d(np.loadtxt(f))
    tt = np.abs(dd[:, C['T_cut']]) / 1e6
    print('  %-10d %-9s %-4s %4d %10.4f %10.4f %10.4f %10.4f'
          % (eta, cond, 'yes' if overturned else 'no', dd.shape[0],
             np.median(tt) / B_NORM, np.median(tt),
             np.median(np.abs(0.5 * dd[:, C['H']] * dd[:, C['K']])),
             np.median(np.abs(dd[:, C['KQn']])) / 1e6))
print()
