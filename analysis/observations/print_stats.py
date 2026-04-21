#!/usr/bin/env python3
"""
Print DP and scaling statistics for reference and sensitivity viscosities.

For each viscosity, reports stats for ALL segments and for scaling < 5 MPa segments:
  - % of segments with scaling < 5 MPa and < 10 MPa
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
    dp, scaling = d[:, 0], d[:, 1]
    valid = ~np.isnan(scaling)
    n_tot = valid.sum()
    n5  = np.sum(np.abs(scaling[valid]) < 5)
    n10 = np.sum(np.abs(scaling[valid]) < 10)
    mask5 = np.abs(scaling) < 5
    dp5 = dp[mask5]
    return {
        'pct5':   n5  / n_tot * 100,
        'pct10':  n10 / n_tot * 100,
        'all_min':  dp.min(),
        'all_mean': dp.mean(),
        'all_max':  dp.max(),
        'dp_min':   dp5.min(),
        'dp_mean':  dp5.mean(),
        'dp_max':   dp5.max(),
        'n_tot': int(n_tot),
        'n5':    int(n5),
    }

hdr = f"{'eta (Pa s)':<14} {'subset':<12} {'<5 MPa':>8} {'<10 MPa':>9}  {'DP min':>7} {'DP mean':>8} {'DP max':>7}   N"
print(hdr)
print('-' * len(hdr))
for v in VISCOSITIES:
    s = stats(v)
    print(f"{v:<14.2e} {'all':<12} {s['pct5']:>7.1f}% {s['pct10']:>8.1f}%"
          f"  {s['all_min']:>6.1f}  {s['all_mean']:>7.1f}  {s['all_max']:>6.1f} MPa"
          f"   {s['n_tot']}")
    print(f"{'':14} {'scaling<5 MPa':<12} {'':>8}  {'':>9}"
          f"  {s['dp_min']:>6.1f}  {s['dp_mean']:>7.1f}  {s['dp_max']:>6.1f} MPa"
          f"   {s['n5']}")
    print()
