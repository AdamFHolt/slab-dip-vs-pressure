#!/usr/bin/env python3
"""
Print DP and Lambda statistics for reference and sensitivity viscosities.

For each viscosity, reports stats for ALL segments and for Lambda < 0.1 segments:
  - % of segments with Lambda < 0.1 and < 0.2
  - DP min, mean, max

Usage:
    python3 print_stats.py
"""

import numpy as np

ALPHA       = 3.28e-5
TM          = 1333.0
DIFFUSIVITY = 8.044e-7
PLATE_THICK = 88000.0
CRUST_THICK = 7000.0
COOLING     = 'plate-cooling'

VISCOSITIES = [2e22, 4e22, 8e22]

BASE = ('text_files/maps.slab{v}.alpha{a}.T{T}.k{k}'
        '.platethick{pt}.crustthick{ct}.{cm}.txt')

def fname(visc):
    return BASE.format(
        v=visc, a=ALPHA, T=TM, k=DIFFUSIVITY,
        pt=PLATE_THICK, ct=CRUST_THICK, cm=COOLING)

def load(visc):
    return np.genfromtxt(fname(visc), delimiter=',', dtype=float)

def stats(visc):
    d = load(visc)
    if d.ndim == 1:
        d = d.reshape(1, -1)
    dp, lam = d[:, 0], d[:, 2]
    valid = ~np.isnan(lam)
    n_tot = valid.sum()
    n1  = np.sum(np.abs(lam[valid]) < 0.1)
    n2  = np.sum(np.abs(lam[valid]) < 0.2)
    mask1 = np.abs(lam) < 0.1
    dp1 = dp[mask1]
    return {
        'pct1':   n1 / n_tot * 100,
        'pct2':   n2 / n_tot * 100,
        'all_min':  dp.min(),
        'all_mean': dp.mean(),
        'all_max':  dp.max(),
        'dp_min':   dp1.min(),
        'dp_mean':  dp1.mean(),
        'dp_max':   dp1.max(),
        'n_tot': int(n_tot),
        'n1':    int(n1),
    }

hdr = f"{'eta (Pa s)':<14} {'subset':<12} {'L<0.1':>8} {'L<0.2':>9}  {'DP min':>7} {'DP mean':>8} {'DP max':>7}   N"
print(hdr)
print('-' * len(hdr))
for v in VISCOSITIES:
    s = stats(v)
    print(f"{v:<14.2e} {'all':<12} {s['pct1']:>7.1f}% {s['pct2']:>8.1f}%"
          f"  {s['all_min']:>6.1f}  {s['all_mean']:>7.1f}  {s['all_max']:>6.1f} MPa"
          f"   {s['n_tot']}")
    print(f"{'':14} {'Lambda<0.1':<12} {'':>8}  {'':>9}"
          f"  {s['dp_min']:>6.1f}  {s['dp_mean']:>7.1f}  {s['dp_max']:>6.1f} MPa"
          f"   {s['n1']}")
    print()
