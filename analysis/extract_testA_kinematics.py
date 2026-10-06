"""
Test A of the shear-vs-bending brief (2026-10-06): does the slab deform in a
through-thickness shear mode (Eq. 6, dQ_s/ds ~ H) or a thin-sheet bending mode
(dQ_s/ds ~ H^3)?

In 2-D, the rotation rate of a material line element with direction e is
e_perp . (e . grad v).  Taking the slab tangent t and normal n_hat at the
midplane (both at the analysis depth, held fixed during the differencing):

    omega_t = n_hat . dv/ds      rotation rate of the midsurface tangent
                                 (= dv_n/ds - K v_s, the bracket kept in Eq. 6)
    omega_n = -dv_s/dn = -a      rotation rate of the material normal
                                 (a = slope of v_s across the slab interior)

Kirchhoff (bending): normals stay normal, omega_n = omega_t, r = -a/omega_t = 1.
Shear mode (paper):  cross-sections translate without rotating, a = 0, r = 0.
The two terms are also the two halves of the shear strain rate,
2 eps_sn = a + omega_t, so tau_sn = eta (a + omega_t) is a check on the
extraction (Test A'); tau_sn / (eta omega_t) is the ratio the brief asks for.

Inputs: csv_outputs/<model>/full.<t>.csv and the extraction text file (dip,
K, H, coverage) from text_files/<src>/, with src = extracted_coverage for the
15 paper models and extracted_varH for the H = 70/130 runs.  Geometry follows
extract_properties.py: midplane point from the horizontal profile at the
analysis depth (max x with c_llith > 0.5), along-s differencing between the
midplane points 10 km shallower and deeper.

Output: text_files/testA/<model>.z<zkm>.txt, one row per timestep, header
inside.  Nothing existing is read for writing or overwritten.

Usage:  python3 extract_testA_kinematics.py <model> [analysis_depth_m=300e3]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

ANALYSIS = os.path.dirname(os.path.abspath(__file__)) + '/'
OUTDIR = ANALYSIS + 'text_files/testA'

YMAX = 1450.e3
DZ_PROF = 1.e3          # half-height of the horizontal-profile band (prof-dz1.0km)
DS_ALONG = 10.e3        # shear-dz10.0: midplane points at z -/+ 10 km
LOCAL_R = 250.e3        # radius of the point cloud used for interpolation
INTERIOR = 0.6          # central fraction of the thickness used for the fits
FACE_CUT = 10.e3        # cut for the Q_s integral, as in extract_properties.py
S_PER_YR = 365.25 * 24 * 3600.
FIRST_TIME = 8
T_START = 11            # first timestep the figures use

SUFFIX = '.z%s.shear-dz10.0.ds10.0.prof-dz1.0km.txt'
VARH = 'H70km', 'H130km'


def horiz_center(x, y, llith, ulith, depth):
    """Midplane, lower and upper face x at one depth, as get_slablocation_from_horiz_prof."""
    yloc = YMAX - depth
    m = (y < yloc + DZ_PROF) & (y > yloc - DZ_PROF)
    ll = m & (llith > 0.5)
    ul = m & (ulith > 0.5)
    if not ll.any() or not ul.any():
        return None
    return x[ll].max(), x[ll].min(), x[ul].max()


def profile_fit(n, vs, lo, hi):
    """Slope of v_s across the interior, (1/s), and the number of points used."""
    half = 0.5 * INTERIOR * (hi - lo)
    mid = 0.5 * (lo + hi)
    sel = (n > mid - half) & (n < mid + half) & np.isfinite(vs)
    if sel.sum() < 5:
        return np.nan, np.nan, sel.sum()
    p = np.polyfit(n[sel], vs[sel], 1)
    return p[0], p[1] + p[0] * mid, sel.sum()


def main():
    model = sys.argv[1]
    depth = float(sys.argv[2]) if len(sys.argv) > 2 else 300.e3
    zkm = '%.1f' % (depth / 1.e3)
    src = 'extracted_varH' if model.endswith(VARH) else 'extracted_coverage'
    txt = np.loadtxt(ANALYSIS + 'text_files/%s/%s' % (src, model) + SUFFIX % zkm)
    os.makedirs(OUTDIR, exist_ok=True)

    rows = []
    for irow in range(T_START - FIRST_TIME, txt.shape[0]):
        time = irow + FIRST_TIME
        dip_deg, K, H_txt, cov = txt[irow, 5], txt[irow, 11], txt[irow, 9], txt[irow, 38]
        csv = ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (model, time)
        if not os.path.exists(csv):
            print('%s t=%d: MISSING CSV' % (model, time), flush=True)
            continue
        d = pd.read_csv(csv, usecols=['velocity:0', 'velocity:1', 'shear_stress:0',
                                      'shear_stress:1', 'shear_stress:4', 'viscosity',
                                      'crust', 'ulith', 'llith', 'Points:0', 'Points:1'])
        x = d['Points:0'].to_numpy(); y = d['Points:1'].to_numpy()
        llith = d['llith'].to_numpy(); ulith = d['ulith'].to_numpy()

        c = horiz_center(x, y, llith, ulith, depth)
        c_sh = horiz_center(x, y, llith, ulith, depth - DS_ALONG)
        c_dp = horiz_center(x, y, llith, ulith, depth + DS_ALONG)
        if c is None or c_sh is None or c_dp is None:
            print('%s t=%d: no slab at %.0f km' % (model, time, depth / 1e3), flush=True)
            continue
        xc, yc = c[0], YMAX - depth
        p_sh = np.array([c_sh[0], YMAX - depth + DS_ALONG])
        p_dp = np.array([c_dp[0], YMAX - depth - DS_ALONG])

        # slab frame at the analysis depth: t down-dip, n_hat = t rotated +90 deg
        th = np.deg2rad(dip_deg)
        t = np.array([np.cos(th), -np.sin(th)])
        nh = np.array([np.sin(th), np.cos(th)])

        # one triangulation of the local cloud, reused for every field
        loc = (x - xc) ** 2 + (y - yc) ** 2 < LOCAL_R ** 2
        vals = np.column_stack([d['velocity:0'].to_numpy()[loc] / S_PER_YR,
                                d['velocity:1'].to_numpy()[loc] / S_PER_YR,
                                d['shear_stress:0'].to_numpy()[loc],
                                d['shear_stress:1'].to_numpy()[loc],
                                d['shear_stress:4'].to_numpy()[loc],
                                d['viscosity'].to_numpy()[loc],
                                (llith + ulith)[loc],
                                d['crust'].to_numpy()[loc]])
        interp = LinearNDInterpolator(np.column_stack([x[loc], y[loc]]), vals)

        # slab-normal profile through the midplane point
        n = np.arange(-0.9 * H_txt, 0.9 * H_txt, 0.5e3)
        pts = np.array([xc, yc])[None, :] + n[:, None] * nh[None, :]
        f = interp(pts)
        vx, vy, sxx, sxy, syy, eta, lith, crust = f.T
        inslab = np.where(lith > 0.5)[0]
        if len(inslab) < 10:
            print('%s t=%d: thin lith profile' % (model, time), flush=True)
            continue
        # the contiguous lith segment containing n = 0
        i0 = np.argmin(np.abs(n))
        if lith[i0] <= 0.5:
            i0 = inslab[np.argmin(np.abs(inslab - i0))]
        lo_i = i0
        while lo_i > 0 and lith[lo_i - 1] > 0.5:
            lo_i -= 1
        hi_i = i0
        while hi_i < len(n) - 1 and lith[hi_i + 1] > 0.5:
            hi_i += 1
        n_lo, n_hi = n[lo_i], n[hi_i]
        H_n = n_hi - n_lo

        vs = vx * t[0] + vy * t[1]
        vn = vx * nh[0] + vy * nh[1]
        # tau_ij n_j projected on t (shear) and on t (normal along s)
        tau_sn = t[0] * (sxx * nh[0] + sxy * nh[1]) + t[1] * (sxy * nh[0] + syy * nh[1])
        tau_ss = t[0] * (sxx * t[0] + sxy * t[1]) + t[1] * (sxy * t[0] + syy * t[1])

        a, vs0, npts = profile_fit(n, vs, n_lo, n_hi)
        _, vn0, _ = profile_fit(n, vn, n_lo, n_hi)
        tau1, tau_ss0, _ = profile_fit(n, tau_ss, n_lo, n_hi)
        half = 0.5 * INTERIOR * H_n
        mid = 0.5 * (n_lo + n_hi)
        sel = (n > mid - half) & (n < mid + half)
        eta_int = np.nanmedian(eta[sel])
        tau_sn_mean = np.nanmean(tau_sn[sel])
        cut = (n > n_lo + FACE_CUT) & (n < n_hi - FACE_CUT) & np.isfinite(tau_sn)
        Q_s = np.trapz(tau_sn[cut], n[cut]) if cut.sum() > 2 else np.nan
        M = np.trapz((n[cut] - mid) * tau_ss[cut], n[cut]) if cut.sum() > 2 else np.nan

        # tangent rotation rate from the midplane velocities 10 km up and down slab
        v_sh = interp(p_sh[None, :])[0, :2]
        v_dp = interp(p_dp[None, :])[0, :2]
        ds_tot = np.linalg.norm(p_sh - np.array([xc, yc])) + np.linalg.norm(p_dp - np.array([xc, yc]))
        dv_ds = (v_dp - v_sh) / ds_tot
        omega_t = nh @ dv_ds
        dvs_ds = t @ dv_ds

        # slope of v_s at the shallow and deep points too (for d(omega_n)/ds)
        a_sd = []
        for p in (p_sh, p_dp):
            pp = p[None, :] + n[:, None] * nh[None, :]
            ff = interp(pp)
            vss = ff[:, 0] * t[0] + ff[:, 1] * t[1]
            li = ff[:, 6] > 0.5
            if li.sum() > 10:
                nn = n[li]
                a_sd.append(profile_fit(n, vss, nn.min(), nn.max())[0])
            else:
                a_sd.append(np.nan)
        da_ds = (a_sd[1] - a_sd[0]) / ds_tot

        r = -a / omega_t if omega_t != 0 else np.nan
        rows.append([time, dip_deg, K, H_txt, H_n, eta_int, vs0 * S_PER_YR * 100, vn0 * S_PER_YR * 100,
                     a, omega_t, r, dvs_ds, a + omega_t, tau_sn_mean, Q_s,
                     tau_sn_mean / (eta_int * omega_t), tau_sn_mean / (eta_int * (a + omega_t)),
                     tau_ss0, tau1, M, a_sd[0], a_sd[1], da_ds, npts, n_lo, n_hi, cov])
        print('%s t=%2d  dip %5.1f  H_n %5.1f km  a=%9.2e  omega_t=%9.2e  r=%6.2f  '
              'tau_sn/(eta*2eps)=%6.2f  tau_sn/(eta*omega)=%6.2f  N=%d'
              % (model, time, dip_deg, H_n / 1e3, a, omega_t, r,
                 tau_sn_mean / (eta_int * (a + omega_t)), tau_sn_mean / (eta_int * omega_t), npts),
              flush=True)

    header = ('0 time | 1 dip_deg | 2 K_1/m | 3 H_txt_m | 4 H_normal_m | 5 eta_int | 6 vs0_cm/yr | 7 vn0_cm/yr | '
              '8 a=dvs/dn_1/s | 9 omega_t=n.dv/ds_1/s | 10 r=-a/omega_t | 11 dvs/ds_1/s | 12 2eps_sn_kin_1/s | '
              '13 tau_sn_mean_Pa | 14 Q_s_N/m | 15 tau_sn/(eta*omega_t) | 16 tau_sn/(eta*2eps_kin) | '
              '17 tau_ss0_Pa | 18 dtau_ss/dn_Pa/m | 19 M_N | 20 a_shallow | 21 a_deep | 22 da/ds_1/(m s) | '
              '23 n_fit_pts | 24 n_lo_m | 25 n_hi_m | 26 coverage')
    np.savetxt(os.path.join(OUTDIR, '%s.z%s.txt' % (model, zkm)), np.array(rows), header=header)


if __name__ == '__main__':
    main()
