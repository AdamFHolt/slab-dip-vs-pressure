"""
Test C of the shear-vs-bending brief: moment balance across the slab.

Integrating the s-momentum equation d(sigma_ss)/ds + d(sigma_sn)/dn + rho g_s = 0
across the slab with weight (n - n_mid) gives, at leading order in Kn,

    Q_s = dM/ds + C,     M = int (n - n_mid) sigma_ss dn,
                         C = (H/2) [tau_sn(top face) + tau_sn(bottom face)]

(the face-traction couple), with sigma_ss = tau_ss - p_dyn.  The dynamic
pressure is used because the lithostatic part is balanced by rho_ref g_s and
the remaining Delta-rho g_s body force is symmetric about the midplane.
This identity holds in either deformation regime; what discriminates is
the partition of Q_s between the internal moment gradient (bending) and the
couple exerted by the mantle on the faces (shear mode), and whether the
measured dQ_s/ds equals d2M/ds2 + dC/ds.

Also Test D: the bending rate two ways, kinematic kdot = d(omega_n)/ds with
omega_n = -dv_s/dn, and from the stress profile kdot = -(dtau_ss/dn)/(2 eta);
and M against the thin-sheet constitutive estimate eta H^3 kdot / 3.

Profiles are taken at five depths, z + k*10 km (k = -2..2), through the
midplane point found as in extract_properties.py, all in the slab frame at
z.  M(s) is fitted with a quadratic, Q_s(s), C(s) and a(s) with lines.
Stresses are converted to tension-positive: ASPECT's shear_stress is minus
2 eta eps (checked by Test A: tau_sn/(eta 2eps) = -0.94), so tau = -S.

Output: text_files/testC/<model>.z<zkm>.txt (new directory).
Usage:  python3 extract_testC_moment.py <model> [analysis_depth_m=300e3]
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

from extract_testA_kinematics import (ANALYSIS, YMAX, DS_ALONG, LOCAL_R, S_PER_YR, FIRST_TIME,
                                      T_START, SUFFIX, VARH, horiz_center, profile_fit)

OUTDIR = ANALYSIS + 'text_files/testC'
KS = (-2, -1, 0, 1, 2)
CUT = 2.5e3          # face cut for the M and Q_s integrals (the composition spike only)
CUT_Q = 10.e3        # the extraction's own cut, for comparison
OUT = 3.e3           # distance outside the face at which the mantle traction is read


def lith_segment(n, lith):
    i0 = np.argmin(np.abs(n))
    ins = np.where(lith > 0.5)[0]
    if len(ins) < 10:
        return None
    if lith[i0] <= 0.5:
        i0 = ins[np.argmin(np.abs(ins - i0))]
    lo = i0
    while lo > 0 and lith[lo - 1] > 0.5:
        lo -= 1
    hi = i0
    while hi < len(n) - 1 and lith[hi + 1] > 0.5:
        hi += 1
    return n[lo], n[hi]


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
        dqds_txt = -txt[irow, 6]                     # paper sign, Pa/m
        csv = ANALYSIS + 'csv_outputs/%s/full.%d.csv' % (model, time)
        if not os.path.exists(csv):
            print('%s t=%d: MISSING CSV' % (model, time), flush=True)
            continue
        d = pd.read_csv(csv, usecols=['velocity:0', 'velocity:1', 'shear_stress:0', 'shear_stress:1',
                                      'shear_stress:4', 'nonadiabatic_pressure', 'viscosity',
                                      'ulith', 'llith', 'Points:0', 'Points:1'])
        x = d['Points:0'].to_numpy(); y = d['Points:1'].to_numpy()
        llith = d['llith'].to_numpy(); ulith = d['ulith'].to_numpy()

        centers = []
        for k in KS:
            c = horiz_center(x, y, llith, ulith, depth + k * DS_ALONG)
            if c is None:
                break
            centers.append(np.array([c[0], YMAX - depth - k * DS_ALONG]))
        if len(centers) < 5:
            print('%s t=%d: no slab at one of the five depths' % (model, time), flush=True)
            continue
        xc, yc = centers[2]
        th = np.deg2rad(dip_deg)
        t = np.array([np.cos(th), -np.sin(th)])
        nh = np.array([np.sin(th), np.cos(th)])

        loc = (x - xc) ** 2 + (y - yc) ** 2 < LOCAL_R ** 2
        vals = np.column_stack([d['velocity:0'].to_numpy()[loc] / S_PER_YR,
                                d['velocity:1'].to_numpy()[loc] / S_PER_YR,
                                -d['shear_stress:0'].to_numpy()[loc],      # tension positive
                                -d['shear_stress:1'].to_numpy()[loc],
                                -d['shear_stress:4'].to_numpy()[loc],
                                d['nonadiabatic_pressure'].to_numpy()[loc],
                                d['viscosity'].to_numpy()[loc],
                                (llith + ulith)[loc]])
        interp = LinearNDInterpolator(np.column_stack([x[loc], y[loc]]), vals)

        n = np.arange(-1.0 * H_txt, 1.0 * H_txt, 0.5e3)
        s_pts, M_full, M_dev, Q_pts, Q10_pts, C_pts, a_pts, H_pts = [], [], [], [], [], [], [], []
        center_diag = None
        s_cum = 0.0
        for i, c in enumerate(centers):
            if i > 0:
                s_cum += np.linalg.norm(c - centers[i - 1])
            f = interp(c[None, :] + n[:, None] * nh[None, :])
            vx, vy, sxx, sxy, syy, pdyn, eta, lith = f.T
            seg = lith_segment(n, lith)
            if seg is None:
                break
            n_lo, n_hi = seg
            mid = 0.5 * (n_lo + n_hi); Hn = n_hi - n_lo
            vs = vx * t[0] + vy * t[1]
            tau_ss = t[0] * (sxx * t[0] + sxy * t[1]) + t[1] * (sxy * t[0] + syy * t[1])
            tau_sn = t[0] * (sxx * nh[0] + sxy * nh[1]) + t[1] * (sxy * nh[0] + syy * nh[1])
            tau_nn = nh[0] * (sxx * nh[0] + sxy * nh[1]) + nh[1] * (sxy * nh[0] + syy * nh[1])
            sig_ss = tau_ss - pdyn
            sig_nn = tau_nn - pdyn
            cut = (n > n_lo + CUT) & (n < n_hi - CUT) & np.isfinite(sig_ss)
            cut10 = (n > n_lo + CUT_Q) & (n < n_hi - CUT_Q) & np.isfinite(tau_sn)
            if cut.sum() < 5:
                break
            s_pts.append(s_cum)
            M_full.append(np.trapz((n[cut] - mid) * sig_ss[cut], n[cut]))
            M_dev.append(np.trapz((n[cut] - mid) * tau_ss[cut], n[cut]))
            Q_pts.append(np.trapz(tau_sn[cut], n[cut]))
            Q10_pts.append(np.trapz(tau_sn[cut10], n[cut10]) if cut10.sum() > 2 else np.nan)
            top = np.interp(n_hi + OUT, n, tau_sn); bot = np.interp(n_lo - OUT, n, tau_sn)
            C_pts.append(0.5 * Hn * (top + bot))
            a_pts.append(profile_fit(n, vs, n_lo, n_hi)[0])
            H_pts.append(Hn)
            if i == 2:
                tau1, tau_ss0, _ = profile_fit(n, tau_ss, n_lo, n_hi)
                half = 0.3 * Hn
                sel = (n > mid - half) & (n < mid + half)
                eta_int = np.nanmedian(eta[sel])
                # sigma_nn should be continuous across the faces; the jump is a sign/convention check
                jump_top = np.interp(n_hi + OUT, n, sig_nn) - np.interp(n_hi - OUT, n, sig_nn)
                jump_bot = np.interp(n_lo + OUT, n, sig_nn) - np.interp(n_lo - OUT, n, sig_nn)
                sig_ss_range = np.nanmax(sig_ss[cut]) - np.nanmin(sig_ss[cut])
                center_diag = (tau1, tau_ss0, eta_int, 0.5 * (abs(jump_top) + abs(jump_bot)), sig_ss_range, top, bot, Hn)
        if len(s_pts) < 5:
            print('%s t=%d: profile failed at one depth' % (model, time), flush=True)
            continue
        s = np.array(s_pts) - s_pts[2]
        pM = np.polyfit(s, M_full, 2); pMd = np.polyfit(s, M_dev, 2)
        pQ = np.polyfit(s, Q_pts, 1); pC = np.polyfit(s, C_pts, 1); pa = np.polyfit(s, a_pts, 1)
        dM_ds, d2M_ds2 = pM[1], 2 * pM[0]
        dMd_ds, d2Md_ds2 = pMd[1], 2 * pMd[0]
        Q0, dQ_ds = pQ[1], pQ[0]
        C0, dC_ds = pC[1], pC[0]
        kdot_kin = -pa[0]                       # d(omega_n)/ds, omega_n = -a
        tau1, tau_ss0, eta_int, nn_jump, ss_range, top, bot, Hn = center_diag
        kdot_str = -tau1 / (2 * eta_int)
        M_est = eta_int * Hn ** 3 * kdot_kin / 3.0

        rows.append([time, dip_deg, K, H_txt, Hn, eta_int, cov,
                     Q0, Q10_pts[2], dM_ds, C0, dM_ds + C0,
                     dQ_ds, dqds_txt, d2M_ds2, dC_ds, d2M_ds2 + dC_ds,
                     M_full[2], M_dev[2], dMd_ds, d2Md_ds2,
                     kdot_kin, kdot_str, M_est, top, bot, nn_jump, ss_range])
        print('%s t=%2d  Q=%7.2f  dM/ds=%7.2f  C=%7.2f  (sum %7.2f) MPa.m/m->MPa | '
              'dQ/ds: fit %6.2f txt %6.2f  d2M+dC %6.2f | kdot kin/str %.2f | M/M_est %.2f | nn_jump/ss_range %.2f'
              % (model, time, Q0 / 1e6 / Hn * Hn / 1e3 * 1e3 / 1e3 if False else Q0 / 1e9, dM_ds / 1e9, C0 / 1e9, (dM_ds + C0) / 1e9,
                 dQ_ds / 1e6, dqds_txt / 1e6, (d2M_ds2 + dC_ds) / 1e6,
                 kdot_kin / kdot_str if kdot_str else np.nan, M_full[2] / M_est if M_est else np.nan,
                 nn_jump / ss_range if ss_range else np.nan), flush=True)

    header = ('0 time | 1 dip | 2 K | 3 H_txt | 4 H_n | 5 eta_int | 6 cov | '
              '7 Q_s(2.5km cut) N/m | 8 Q_s(10km cut) | 9 dM/ds | 10 C couple | 11 dM/ds+C | '
              '12 dQ/ds fit Pa/m | 13 dQ/ds txt (paper sign) | 14 d2M/ds2 | 15 dC/ds | 16 d2M/ds2+dC/ds | '
              '17 M_full N | 18 M_dev | 19 dM_dev/ds | 20 d2M_dev/ds2 | '
              '21 kdot_kin 1/(m s) | 22 kdot_stress | 23 M_est=eta H^3 kdot_kin/3 | '
              '24 tau_sn top face | 25 tau_sn bottom face | 26 sigma_nn face jump | 27 sigma_ss range')
    np.savetxt(os.path.join(OUTDIR, '%s.z%s.txt' % (model, zkm)), np.array(rows), header=header)


if __name__ == '__main__':
    main()
