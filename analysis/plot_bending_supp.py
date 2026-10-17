"""
Supplement figure: the slab deforms as a bending sheet, so dQ_s/ds comes from
the internal bending moment (thin-sheet route, dQ_s/ds ~ eta H^3 K v_c /
(3 L_K^3)) and not from through-thickness shear (Eq. 6 of the submitted
version, ~ eta H).

  a) tau_sn(n) across the slab for one model and output, against the
     shear-mode estimate eta*omega (Eq. 6 of the submitted version)
  b) cross-section rotation rate -dv_s/dn against the midsurface rotation
     rate omega for every output (Test A): bending on 1:1, shear mode on 0
  c) Q_s split into the moment gradient dM/ds and the face-traction couple C (Test C)
  d) |M| against the thin-sheet moment |eta h^3 kdot / 3| (Test C)
(2x2 since 2026-10-17: the v_s(n) and tau_ss(n) profiles were dropped, as
the one-output versions of b) and d))

Model set and selection as Fig. 5 (plot_forces.varH.py ... alleta nodrho71):
five strengths at H = 100, three at H = 70/130, every plate configuration,
t >= 11, covered rows only, the same overturned set (hollow, left out of the
statistics) and the fixed-SP eta' = 50 trim.  Inputs: text_files/testA and
text_files/testC (extract_testA_kinematics.py, extract_testC_moment.py).

Output: plots/DP-comparisons/compilations/manuscript_figs/bending_supp.z<zkm>[.<model>.t<t>].png/pdf
Usage:  python3 plot_bending_supp.py [depth_m=300e3] [model:time]
        (default profile: the reference model new2 at output 20)
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import LinearNDInterpolator

import extract_testA_kinematics as A

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams['font.size'] = 7
mpl.rcParams.update({'axes.labelsize': 9, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5})

args = [a for a in sys.argv[1:] if ':' not in a]
depth = float(args[0]) if args else 300.e3
zkm = '%.1f' % (depth / 1.e3)
prof = [a for a in sys.argv[1:] if ':' in a]
PROF_MODEL, PROF_T = (prof[0].split(':')[0], int(prof[0].split(':')[1])) if prof \
    else ('2D_compositional_subd_lower-res_new2', 20)

# --- model set, as plot_forces.varH.py
H100 = {'free':    {50: 'lower-res_new_50plates', 250: 'lower-res_new_250plates', 375: 'lower-res_new_375plates',
                    500: 'lower-res_new2', 1000: 'lower-res_new_1000plates'},
        'FixedSP': {50: 'lower-res_new_FixedSP_50plates', 250: 'lower-res_new_FixedSP_250plates',
                    375: 'lower-res_new_FixedSP_375plates', 500: 'FixedSP_lower-res_new2',
                    1000: 'lower-res_new_FixedSP_1000plates2'},
        'FixedOP': {50: 'lower-res_new_FixedOP_50plates', 250: 'lower-res_new_FixedOP_250plates',
                    375: 'lower-res_new_FixedOP_375plates', 500: 'FixedOP_lower-res_new',
                    1000: 'lower-res_new_FixedOP_1000plates'}}


def name(bc, eta, H):
    if H == 100:
        return '2D_compositional_subd_' + H100[bc][eta]
    return '2D_compositional_subd_lower-res_new_%s%dplates_H%dkm' % ('' if bc == 'free' else bc + '_', eta, H)


BCS, ETAS, HS = ['free', 'FixedSP', 'FixedOP'], [50, 250, 375, 500, 1000], [70, 100, 130]
colors = {50: 'tan', 250: 'peru', 375: 'firebrick', 500: 'maroon', 1000: 'black'}
bc_marker = {'FixedSP': 'v', 'free': 'o', 'FixedOP': '^'}
bc_label = {'FixedSP': 'fixed SP', 'free': 'free', 'FixedOP': 'fixed OP'}
SIZE = {70: 5, 100: 9, 130: 15}                 # marker area by H, as Fig. 5's insets
OVERTURNED = {('free', 1000, 100), ('free', 500, 130), ('free', 1000, 130),
              ('FixedOP', 500, 100), ('FixedOP', 1000, 100),
              ('FixedOP', 250, 130), ('FixedOP', 500, 130), ('FixedOP', 1000, 130),
              ('FixedOP', 1000, 70), ('FixedOP', 375, 100)}
LAST_OUTPUT = {('FixedSP', 50, 100): 22}


def load(test, bc, eta, H, cov_col):
    f = 'text_files/test%s/%s.z%s.txt' % (test, name(bc, eta, H), zkm)
    if not os.path.exists(f):
        return None
    d = np.loadtxt(f, ndmin=2)
    d = d[(d[:, 0] >= 11) & (d[:, cov_col] == 1)]
    return d[d[:, 0] <= LAST_OUTPUT.get((bc, eta, H), 1e9)]


keys = [(bc, e, H) for bc in BCS for e in ETAS for H in HS if H == 100 or e in (250, 500, 1000)]
TA = {k: load('A', *k, 26) for k in keys}
TC = {k: load('C', *k, 6) for k in keys}
missing = [k for k in keys if TA[k] is None or TC[k] is None]
if missing:
    print('missing Test A/C files:', missing)
keys = [k for k in keys if k not in missing]
normal = [k for k in keys if k not in OVERTURNED]

CUT = 2.5e3      # face cut of the Test C integrals (extract_testC_moment.py)


def M_thin(d):
    """Thin-sheet moment eta h^3 kdot / 3 over the same trimmed thickness h = H_n - 2 CUT
    as the measured M, with kdot = d(dv_s/dn)/ds (Test C's col 21 is minus this).  With the
    full stress in plane strain (sigma_nn ~ 0 across the slab) sigma_ss = 4 eta n kdot, so
    M = int n sigma_ss dn = eta h^3 kdot / 3, tension positive."""
    return -d[:, 23] * ((d[:, 4] - 2 * CUT) / d[:, 4]) ** 3


# --- statistics for the caption (normal runs)
r_all = np.concatenate([TA[k][:, 10] for k in normal])
r_all = r_all[np.isfinite(r_all)]
print('Test A, %s km: %d normal outputs from %d models; r median %.2f, IQR %.2f-%.2f, '
      '%.0f %% > 0.7, %.0f %% < 0.3' % (zkm, len(r_all), len(normal), np.median(r_all),
      *np.percentile(r_all, [25, 75]), 100 * np.mean(r_all > 0.7), 100 * np.mean(r_all < 0.3)))
for H in HS:
    rH = np.concatenate([TA[k][:, 10] for k in normal if k[2] == H])
    print('  H = %3d: median r %.2f (N = %d)' % (H, np.nanmedian(rH), np.isfinite(rH).sum()))
for bc in BCS:
    rB = np.concatenate([TA[k][:, 10] for k in normal if k[0] == bc])
    print('  %-7s: median r %.2f (N = %d)' % (bc, np.nanmedian(rB), np.isfinite(rB).sum()))
tA = np.concatenate([TA[k][:, 15] for k in normal])
print('  outside -0.5 < r < 2 (piled into the end bins of d): %d of %d'
      % (np.sum((r_all <= -0.5) | (r_all >= 2)), len(r_all)))
wN = np.concatenate([TA[k][:, 9] for k in normal]); aN = -np.concatenate([TA[k][:, 8] for k in normal])
print('  -dv_s/dn on omega: slope through the origin %.2f, r = %.3f; outside panel d: %d'
      % (np.sum(wN * aN) / np.sum(wN ** 2), np.corrcoef(wN, aN)[0, 1],
         np.sum((wN * 1e6 * A.S_PER_YR < -0.1) | (wN * 1e6 * A.S_PER_YR > 0.03))))
print('  tau_sn / (eta omega_t): median %.2f' % np.nanmedian(tA))

C = np.vstack([TC[k] for k in normal])
big = np.abs(C[:, 7]) > 1e11                      # |Q_s| > 100 GN/m, so the ratios are defined
print('Test C, %s km: %d normal outputs (%d with |Q_s| > 100 GN/m)' % (zkm, len(C), big.sum()))
print('  dM/ds / Q_s median %.2f (IQR %.2f-%.2f);  C / Q_s median %.2f;  (dM/ds + C)/Q_s median %.2f'
      % (np.median(C[big, 9] / C[big, 7]), *np.percentile(C[big, 9] / C[big, 7], [25, 75]),
         np.median(C[big, 10] / C[big, 7]), np.median(C[big, 11] / C[big, 7])))
pM = np.sum(C[:, 7] * C[:, 9]) / np.sum(C[:, 7] ** 2)
pC = np.sum(C[:, 7] * C[:, 10]) / np.sum(C[:, 7] ** 2)
print('  slopes through the origin on Q_s: dM/ds %.2f, C %.2f' % (pM, pC))
Mr = C[:, 17] / M_thin(C)
print('  M / (eta h^3 kdot/3), h = H_n - 2 x 2.5 km: median %.2f (IQR %.2f-%.2f); eta\' >= 250 only: %.2f'
      % (np.nanmedian(Mr), *np.nanpercentile(Mr, [25, 75]),
         np.nanmedian(np.concatenate([TC[k][:, 17] / M_thin(TC[k]) for k in normal if k[1] >= 250]))))
print('  M and eta h^3 kdot/3 same sign in %.1f %% of outputs; M < 0 in %.0f %%'
      % (100 * np.mean(np.sign(C[:, 17]) == np.sign(M_thin(C))), 100 * np.mean(C[:, 17] < 0)))
print('  same with the full H_n: median %.2f' % np.nanmedian(C[:, 17] / (-C[:, 23])))
for H in HS:
    MH = np.concatenate([np.abs(TC[k][:, 17]) for k in normal if k[2] == H])
    print('  H = %3d: median |M| %.2e N' % (H, np.median(MH)))
cut_fac = {H: ((H * 0.87 - 2 * 2.5) / (H * 0.87)) ** 3 for H in HS}   # rough, measured H ~ 0.87 nominal
print('  (the 2.5 km face cuts alone remove ~%s of a linear-stress moment)'
      % ', '.join('%.0f %%' % (100 * (1 - cut_fac[H])) for H in HS))

# --- profile through the slab (as plot_testA_profiles.py, tension positive)
src = 'extracted_varH' if '_H' in PROF_MODEL else 'extracted_coverage'
txt = np.loadtxt(A.ANALYSIS + 'text_files/%s/%s' % (src, PROF_MODEL) + A.SUFFIX % zkm)
dip_deg, H_txt = txt[PROF_T - A.FIRST_TIME, 5], txt[PROF_T - A.FIRST_TIME, 9]
d = pd.read_csv(A.ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (PROF_MODEL, PROF_T),
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
                        -d['shear_stress:0'].to_numpy()[loc], -d['shear_stress:1'].to_numpy()[loc],
                        -d['shear_stress:4'].to_numpy()[loc], d['viscosity'].to_numpy()[loc], (llith + ulith)[loc]])
interp = LinearNDInterpolator(np.column_stack([x[loc], y[loc]]), vals)
n = np.arange(-1.0 * H_txt, 1.0 * H_txt, 0.5e3)
vx, vy, sxx, sxy, syy, eta, lith = interp(np.array([xc, yc])[None, :] + n[:, None] * nh[None, :]).T
vs = vx * t[0] + vy * t[1]
tau_sn = t[0] * (sxx * nh[0] + sxy * nh[1]) + t[1] * (sxy * nh[0] + syy * nh[1])
tau_ss = t[0] * (sxx * t[0] + sxy * t[1]) + t[1] * (sxy * t[0] + syy * t[1])
p_sh = np.array([c_sh[0], yc + A.DS_ALONG]); p_dp = np.array([c_dp[0], yc - A.DS_ALONG])
ds_tot = np.linalg.norm(p_sh - [xc, yc]) + np.linalg.norm(p_dp - [xc, yc])
omega_t = nh @ ((interp(p_dp[None, :])[0, :2] - interp(p_sh[None, :])[0, :2]) / ds_tot)
ins = np.where(lith > 0.5)[0]
n_lo, n_hi = n[ins.min()], n[ins.max()]
mid = 0.5 * (n_lo + n_hi)
a, vs0, _ = A.profile_fit(n, vs, n_lo, n_hi)
half = 0.5 * A.INTERIOR * (n_hi - n_lo)
eta_int = np.nanmedian(eta[(n > mid - half) & (n < mid + half)])
nk = (n - mid) / 1e3
cm = 100 * A.S_PER_YR
print('profile %s t=%d: dip %.0f, H_n %.0f km, r = %.2f, eta*omega_t = %.0f MPa, interior tau_sn %.1f MPa'
      % (PROF_MODEL, PROF_T, dip_deg, (n_hi - n_lo) / 1e3, -a / omega_t, eta_int * omega_t / 1e6,
         np.nanmean(tau_sn[(n > mid - half) & (n < mid + half)]) / 1e6))

# --- figure: a) tau_sn profile, b) rotation rates, c) Q_s partition, d) moment
fig, axs = plt.subplots(2, 2, figsize=(6.3, 6.3))
(axS, axR), (axQ, axM) = axs
inside = (n >= n_lo) & (n <= n_hi)
# a) three stacked profiles in the top-left slot, sharing n: v_s (-> b), tau_sn
# (the shear-mode claim) and tau_ss (the moment of c, d)
slot = axS.get_subplotspec()
axS.remove()
sub = slot.subgridspec(3, 1, hspace=0.12)
axA = [fig.add_subplot(sub[0])]
axA += [fig.add_subplot(sub[i], sharex=axA[0]) for i in (1, 2)]
axV, axS, axT = axA
for ax in axA:
    ax.axvspan((n_lo - mid) / 1e3, (n_hi - mid) / 1e3, color='0.93', zorder=0, linewidth=0)
    ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
    ax.set_xlim(nk[0], nk[-1])
for ax in (axV, axS):
    ax.tick_params(labelbottom=False)
axT.set_xlabel(r'$n$   [km]')
PL = dict(fontsize=5.5, frameon=False, handlelength=1.6, borderpad=0.2, labelspacing=0.2, handletextpad=0.4)
axV.plot(nk, (vs0 - (n - mid) * omega_t) * cm, color='tab:red', linestyle='--', linewidth=0.9, zorder=2,
         label=r'bending: $v_{s0} - n\,\omega$')
axV.plot(nk, np.full_like(n, vs0) * cm, color='tab:blue', linestyle=':', linewidth=1.1, zorder=2,
         label=r'shear mode: $v_{s0}$')
axV.plot(nk, vs * cm, color='black', linewidth=1.1, zorder=3)
axV.set_ylabel(r'$v_s$ [cm/yr]', fontsize=7.5)
# y-range set by the slab interior, so the tilt shows; the mantle runs off the axes
v_in = np.concatenate([vs[inside], vs0 - (n[inside] - mid) * omega_t]) * cm
pad = 0.35 * (v_in.max() - v_in.min())
axV.set_ylim(v_in.min() - pad - 0.1, v_in.max() + pad + 0.1)
axV.legend(loc='upper left', bbox_to_anchor=(0.09, 1.0), **PL)
est = np.where(inside, eta_int * omega_t, 0.)
axS.plot(nk, est / 1e6, color='tab:blue', linestyle=':', linewidth=1.1, zorder=2, label=r'shear mode: $\eta\,\omega$')
axS.plot(nk, tau_sn / 1e6, color='black', linewidth=1.1, zorder=3)
axS.axhline(0, color='lightgray', linewidth=1, zorder=0)
lim = 1.3 * abs(eta_int * omega_t) / 1e6
axS.set_ylim(-lim, lim)
axS.set_ylabel(r'$\tau_{sn}$ [MPa]', fontsize=7.5)
axS.legend(loc='lower center' if omega_t > 0 else 'upper center', **PL)
axT.plot(nk, tau_ss / 1e6, color='black', linewidth=1.1, zorder=3)
axT.axhline(0, color='lightgray', linewidth=1, zorder=0)
axT.set_ylabel(r'$\tau_{ss}$ [MPa]', fontsize=7.5)
from matplotlib.ticker import MaxNLocator
axT.set_xticks([-50, 0, 50])
for ax in axA:
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_label_coords(-0.16, 0.5)

# d) cross-section rotation rate -dv_s/dn against the midsurface rotation rate omega,
# every output (a scatter rather than a histogram of r = their ratio, which blows up
# where omega passes through zero): bending on the 1:1 line, shear mode on y = 0
MYR = 1e6 * A.S_PER_YR
WL = (-0.1, 0.03)
axR.axhline(0, color='tab:blue', linestyle=':', linewidth=1.2, zorder=1)
axR.text(-0.095, 0.001, 'shear mode', ha='left', va='bottom', fontsize=6.5, color='tab:blue')
axR.plot(WL, WL, color='tab:red', linestyle='--', linewidth=1, zorder=1)
axR.text(-0.068, -0.062, 'bending', ha='center', va='bottom', fontsize=6.5, color='tab:red',
         rotation=45, rotation_mode='anchor')
axR.axvline(0, color='lightgray', linewidth=1, zorder=0)
for k in keys:
    bc, e, H = k
    w, ma = TA[k][:, 9] * MYR, -TA[k][:, 8] * MYR
    if k in OVERTURNED:
        axR.scatter(w, ma, s=SIZE[H], facecolors='none', edgecolor=colors[e], marker=bc_marker[bc],
                    linewidth=0.4, alpha=0.7, zorder=2)
    else:
        axR.scatter(w, ma, s=SIZE[H], color=colors[e], edgecolor='black', marker=bc_marker[bc],
                    linewidth=0.25, zorder=3)
axR.set_xlim(WL); axR.set_ylim(WL)
axR.set_aspect('equal', adjustable='box')
axR.set_xlabel(r'$\omega = \partial v_n/\partial s - K v_s$   [Myr$^{-1}$]')
axR.set_ylabel(r'$-\partial v_s/\partial n$   [Myr$^{-1}$]')

# e) partition of Q_s, f) moment against the thin-sheet estimate
QL = (-400, 2500)
ML = (1e15, 2e18)
for k in keys:
    bc, e, H = k
    dC = TC[k]
    m = bc_marker[bc]
    Q, dM, Cc = dC[:, 7] / 1e9, dC[:, 9] / 1e9, dC[:, 10] / 1e9
    Mm, Me = np.abs(dC[:, 17]), np.abs(M_thin(dC))      # same sign in nearly every output (printed)
    if k in OVERTURNED:
        kw = dict(facecolors='none', edgecolor=colors[e], marker=m, linewidth=0.4, alpha=0.7, zorder=2)
    else:
        kw = dict(color=colors[e], edgecolor='black', marker=m, linewidth=0.25, zorder=3)
    axQ.scatter(Q, dM, s=SIZE[H], **kw)
    axQ.scatter(Q, Cc, s=SIZE[H] * 0.8, marker='x', color=colors[e], linewidth=0.5,
                alpha=0.4 if k in OVERTURNED else 0.8, zorder=2)
    axM.scatter(Me, Mm, s=SIZE[H], **kw)
axM.set_xscale('log'); axM.set_yscale('log')
for ax, (lo, hi), xl, yl in ((axQ, QL, r'$Q_s$   [GN/m]', r'$dM/ds$,  $C$   [GN/m]'),
                              (axM, ML, r'$|\eta H^3 \dot\kappa / 3|$   [N]', r'$|M|$   [N]')):
    ax.plot([lo, hi], [lo, hi], color='gray', linestyle='--', linewidth=1, zorder=1)
    if ax is axQ:
        ax.axhline(0, color='lightgray', linewidth=1, zorder=0)
        ax.axvline(0, color='lightgray', linewidth=1, zorder=0)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_aspect('equal', adjustable='box')
# key for c): the two parts of Q_s (colours, markers and sizes as the strip)
axQ.legend(handles=[Line2D([], [], marker='o', color='gray', markeredgecolor='black', markeredgewidth=0.25,
                           linestyle='', markersize=4, label=r'$dM/ds$'),
                    Line2D([], [], marker='x', color='gray', markeredgewidth=0.6, linestyle='', markersize=4,
                           label=r'$C$')],
           loc='upper left', bbox_to_anchor=(0.01, 0.90), fontsize=6.5, frameon=True, fancybox=False,
           edgecolor='black', facecolor='white', framealpha=1, borderpad=0.5, handletextpad=0.3).get_frame().set_linewidth(0.5)
for ax in list(axA) + [axR, axQ, axM]:
    ax.grid(True, color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
for ax, lab in zip([axV, axR, axQ, axM], 'abcd'):
    ax.text(0.04, 0.96, lab + ')', transform=ax.transAxes, ha='left', va='top', fontsize=9, zorder=10,
            bbox=dict(facecolor='white', edgecolor='none', pad=1))

# legend strip for e), f), as Fig. 5
header = lambda s: Line2D([], [], linestyle='', marker='', label=s)
hs = [header("viscosity ratio η':")] + \
     [Line2D([], [], marker='o', color=colors[e], markeredgecolor='black', markeredgewidth=0.25, linestyle='',
             markersize=4, label='%d' % e) for e in ETAS] + [header('plates:')] + \
     [Line2D([], [], marker=bc_marker[b], color='gray', markeredgecolor='black', markeredgewidth=0.25,
             linestyle='', markersize=4, label=bc_label[b]) for b in ['FixedSP', 'free', 'FixedOP']] + \
     [Line2D([], [], marker='o', markerfacecolor='none', markeredgecolor='gray', linestyle='', markersize=4,
             label='overturned'), header('thickness:')] + \
     [Line2D([], [], marker='o', color='gray', markeredgecolor='black', markeredgewidth=0.25, linestyle='',
             markersize=np.sqrt(SIZE[h]) * 1.0, label='%d km' % h) for h in HS]
plt.tight_layout(h_pad=1.5, w_pad=1.5)
fig.subplots_adjust(top=0.91)
leg = fig.legend(handles=hs, ncol=len(hs), loc='lower center', bbox_to_anchor=(0.5, 0.93), fontsize=6.5,
                 frameon=True, fancybox=False, edgecolor='black', facecolor='white', framealpha=1,
                 borderpad=0.5, handletextpad=0.1, handlelength=1.0, columnspacing=0.5)
leg.get_frame().set_linewidth(0.5)
for h_, tx in zip(hs, leg.get_texts()):
    if h_.get_marker() in ('', 'None', None):
        tx.set_fontweight('bold')

out_dir = 'plots/DP-comparisons/compilations/manuscript_figs/'
tag = '' if not prof else '.%s.t%d' % (PROF_MODEL.replace('2D_compositional_subd_', ''), PROF_T)
out = out_dir + 'bending_supp.z%s%s' % (zkm, tag)
plt.savefig(out + '.png', dpi=300, bbox_inches='tight')
plt.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
