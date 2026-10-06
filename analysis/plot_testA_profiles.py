"""
Test B of the shear-vs-bending brief: slab-normal profiles of v_s, v_n,
tau_sn and tau_ss at the analysis depth for a few (model, timestep) pairs,
with the Kirchhoff prediction v_s(n) = v_s0 - n*omega_t overlaid.

Uniform v_s(n) + uniform tau_sn -> shear mode.  Linear v_s(n) with slope
-omega_t, small tau_sn and linear tau_ss -> bending.

Figure: plots/tmp/testA_profiles.z<zkm>.png (new file).
Usage:  python3 plot_testA_profiles.py <model:time> [<model:time> ...] [--depth 300e3]
"""
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator

import extract_testA_kinematics as A

args = [a for a in sys.argv[1:] if ':' in a]
depth = float(sys.argv[sys.argv.index('--depth') + 1]) if '--depth' in sys.argv else 300.e3
zkm = '%.1f' % (depth / 1.e3)

fig, axs = plt.subplots(len(args), 4, figsize=(13, 2.8 * len(args)), squeeze=False)
for row, arg in enumerate(args):
    model, time = arg.split(':'); time = int(time)
    src = 'extracted_varH' if model.endswith(A.VARH) else 'extracted_coverage'
    txt = np.loadtxt(A.ANALYSIS + 'text_files/%s/%s' % (src, model) + A.SUFFIX % zkm)
    dip_deg, K, H_txt = txt[time - A.FIRST_TIME, 5], txt[time - A.FIRST_TIME, 11], txt[time - A.FIRST_TIME, 9]
    d = pd.read_csv(A.ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (model, time),
                    usecols=['velocity:0', 'velocity:1', 'shear_stress:0', 'shear_stress:1',
                             'shear_stress:4', 'viscosity', 'ulith', 'llith', 'Points:0', 'Points:1'])
    x = d['Points:0'].to_numpy(); y = d['Points:1'].to_numpy()
    llith = d['llith'].to_numpy(); ulith = d['ulith'].to_numpy()
    c = A.horiz_center(x, y, llith, ulith, depth)
    c_sh = A.horiz_center(x, y, llith, ulith, depth - A.DS_ALONG)
    c_dp = A.horiz_center(x, y, llith, ulith, depth + A.DS_ALONG)
    xc, yc = c[0], A.YMAX - depth
    th = np.deg2rad(dip_deg)
    t = np.array([np.cos(th), -np.sin(th)]); nh = np.array([np.sin(th), np.cos(th)])
    loc = (x - xc) ** 2 + (y - yc) ** 2 < A.LOCAL_R ** 2
    vals = np.column_stack([d['velocity:0'].to_numpy()[loc] / A.S_PER_YR, d['velocity:1'].to_numpy()[loc] / A.S_PER_YR,
                            d['shear_stress:0'].to_numpy()[loc], d['shear_stress:1'].to_numpy()[loc],
                            d['shear_stress:4'].to_numpy()[loc], d['viscosity'].to_numpy()[loc], (llith + ulith)[loc]])
    interp = LinearNDInterpolator(np.column_stack([x[loc], y[loc]]), vals)
    n = np.arange(-1.2 * H_txt, 1.2 * H_txt, 0.5e3)
    f = interp(np.array([xc, yc])[None, :] + n[:, None] * nh[None, :])
    vx, vy, sxx, sxy, syy, eta, lith = f.T
    vs = vx * t[0] + vy * t[1]; vn = vx * nh[0] + vy * nh[1]
    tau_sn = t[0] * (sxx * nh[0] + sxy * nh[1]) + t[1] * (sxy * nh[0] + syy * nh[1])
    tau_ss = t[0] * (sxx * t[0] + sxy * t[1]) + t[1] * (sxy * t[0] + syy * t[1])
    p_sh = np.array([c_sh[0], yc + A.DS_ALONG]); p_dp = np.array([c_dp[0], yc - A.DS_ALONG])
    ds_tot = np.linalg.norm(p_sh - [xc, yc]) + np.linalg.norm(p_dp - [xc, yc])
    omega_t = nh @ ((interp(p_dp[None, :])[0, :2] - interp(p_sh[None, :])[0, :2]) / ds_tot)
    ins = np.where(lith > 0.5)[0]; n_lo, n_hi = n[ins.min()], n[ins.max()]
    mid = 0.5 * (n_lo + n_hi)
    a, vs0, _ = A.profile_fit(n, vs, n_lo, n_hi)
    cm = 100 * A.S_PER_YR
    lab = model.replace('2D_compositional_subd_', '').replace('lower-res_', '')
    for ax in axs[row]:
        ax.axvspan(n_lo / 1e3, n_hi / 1e3, color='0.92', zorder=0)
        ax.axvline(mid / 1e3, color='0.7', lw=0.5)
    axs[row, 0].plot(n / 1e3, vs * cm, 'k-', lw=1.2, label=r'$v_s$')
    axs[row, 0].plot(n / 1e3, (vs0 - (n - mid) * omega_t) * cm, 'r--', lw=1, label=r'Kirchhoff: $v_{s0} - n\,\omega_t$')
    axs[row, 0].plot(n / 1e3, np.full_like(n, vs0) * cm, 'b:', lw=1, label=r'shear mode: $v_{s0}$')
    axs[row, 0].set_ylabel('cm/yr'); axs[row, 0].legend(fontsize=6, frameon=False)
    axs[row, 0].set_title('%s  t=%d  dip %.0f  r=%.2f' % (lab, time, dip_deg, -a / omega_t), fontsize=8)
    axs[row, 1].plot(n / 1e3, vn * cm, 'k-', lw=1.2); axs[row, 1].set_title(r'$v_n$', fontsize=8); axs[row, 1].set_ylabel('cm/yr')
    axs[row, 2].plot(n / 1e3, tau_sn / 1e6, 'k-', lw=1.2, label=r'$\tau_{sn}$ (ASPECT sign)')
    axs[row, 2].plot(n / 1e3, -eta * omega_t / 1e6, 'b:', lw=1, label=r'$-\eta\,\omega_t$ (Eq. 6 estimate)')
    axs[row, 2].set_ylabel('MPa'); axs[row, 2].legend(fontsize=6, frameon=False); axs[row, 2].set_title(r'$\tau_{sn}$', fontsize=8)
    axs[row, 3].plot(n / 1e3, tau_ss / 1e6, 'k-', lw=1.2); axs[row, 3].set_title(r'$\tau_{ss}$', fontsize=8); axs[row, 3].set_ylabel('MPa')
    for ax in axs[row]:
        ax.set_xlim(-1.2 * H_txt / 1e3, 1.2 * H_txt / 1e3)
    axs[row, 2].set_ylim(-1.5 * np.nanmax(np.abs(eta * omega_t)) / 1e6, 1.5 * np.nanmax(np.abs(eta * omega_t)) / 1e6)
for ax in axs[-1]:
    ax.set_xlabel('n [km] (lower face -> upper face)')
plt.tight_layout()
out = 'plots/tmp/testA_profiles.z%s.png' % zkm
plt.savefig(out, dpi=130)
print('wrote', out)
