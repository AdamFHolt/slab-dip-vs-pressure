#!/usr/bin/env python3
"""
dQ/ds diagnostic script.

Tests the direct Equation 8 computation (Supporting Text S1):

    dQ/ds ≈ ηH [2K (dvs/ds) + vs (dK/ds)]

against:
  - The measured dQ/ds (col 6, negated to paper sign convention)
  - The paper's empirical Eq. 11:  dQ/ds ≈ 0.1 η K v_c

Produces two figures:
  fig1_formula_comparison  — Eq.8 full, velocity term alone, Eq.11 vs measured
  fig2_calibration         — prefactor vs K, K-filtered histogram, calibrated scatter

Sign convention throughout: paper's dQ/ds > 0 = force resisting slab-normal buoyancy.
  code col 6 (slab_stress_term) = -dQ/ds_paper, so meas_paper = -col6.

Usage (run from analysis/):
    python3 plot_dqds_diagnostic.py [text_files/TESTC/*.txt ...]

Outputs to:  plots/tmp/
"""

import os
import re
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Column indices (extract_properties.py, 30 cols, 0-indexed) ───────────────
COL_TIME    = 0
COL_DP_MOD  = 3    # Pa   — measured across-slab ΔP
COL_DP_ANAL = 4    # Pa   — analytical buoyancy term
COL_DIP     = 5    # deg  — slab dip at analysis depth
COL_DQDS    = 6    # Pa   — slab_stress_term = -dQ/ds_paper
COL_H       = 9    # m    — slab-normal thickness
COL_K       = 11   # rad/m
COL_DK_SG   = 12   # rad/m²  — SG-smoothed dK/ds from full profile
COL_K_SHALL = 13   # rad/m   — K at analysis_depth - dz
COL_K_DEEP  = 14   # rad/m   — K at analysis_depth + dz
COL_VC      = 19   # cm/yr   — convergence velocity
COL_VS      = 21   # m/yr    — slab-parallel velocity at depth
COL_DVDS    = 22   # yr⁻¹    — dvs/ds  (m/yr)/m
COL_LV      = 23   # m       — L_v = |vs / (dvs/ds)|
COL_LK      = 24   # m       — L_K = |K / (dK/ds)|  from SG
COL_ETA     = 26   # Pa·s    — slab viscosity at depth

SPY = 365.25 * 24.0 * 3600.0          # seconds per year
CMYR_TO_MS = 0.01 / SPY               # cm/yr  → m/s


# ── Helpers ──────────────────────────────────────────────────────────────────

def ensure_2d(arr):
    return arr.reshape(1, -1) if arr.ndim == 1 else arr


def parse_shear_dz_km(fpath):
    """Extract analysis_depth_dz in km from filename token  shear-dz<value>."""
    m = re.search(r'shear-dz([\d.]+?)\.(?:ds|km|[a-zA-Z])', os.path.basename(fpath))
    return float(m.group(1)) if m else None


# ── Data loading ─────────────────────────────────────────────────────────────

def load_records(files):
    """
    Return list of per-timestep dicts with all quantities needed for diagnostics.
    """
    records = []
    for fpath in files:
        dz_km = parse_shear_dz_km(fpath)
        if dz_km is None:
            print(f"  skip (no shear-dz in name): {fpath}")
            continue
        dz_m = dz_km * 1e3

        try:
            data = ensure_2d(np.loadtxt(fpath))
        except Exception as e:
            print(f"  skip (load error): {fpath} — {e}")
            continue

        if data.shape[1] < 27:
            print(f"  skip (only {data.shape[1]} cols): {fpath}")
            continue

        model_label = os.path.basename(fpath).replace('.txt', '')
        for prefix in ('2D_compositional_subd_lower-res_',
                       '2D_compositional_subd_'):
            if model_label.startswith(prefix):
                model_label = model_label[len(prefix):]
                break
        model_label = re.sub(r'\.z[\d.]+\.shear.*', '', model_label)

        for row in data:
            dip_deg = row[COL_DIP]
            dip_rad = np.deg2rad(dip_deg)
            if dip_rad <= 0 or not np.isfinite(dip_rad):
                continue

            H        = row[COL_H]
            K        = row[COL_K]
            dK_sg    = row[COL_DK_SG]
            K_shall  = row[COL_K_SHALL]
            K_deep   = row[COL_K_DEEP]
            vs_myr   = row[COL_VS]
            dvds_yr  = row[COL_DVDS]
            eta      = row[COL_ETA]
            vc_cmyr  = row[COL_VC]
            dqds_col = row[COL_DQDS]
            Lv       = row[COL_LV]
            Lk_sg    = row[COL_LK]

            vs_si   = vs_myr  / SPY
            dvds_si = dvds_yr / SPY
            vc_si   = vc_cmyr * CMYR_TO_MS

            meas_pa = -dqds_col   # paper sign: +ve = resists buoyancy

            # 3-point dK/ds from K_shall / K_deep
            ds_total = 2.0 * dz_m / np.sin(dip_rad)
            if (ds_total > 0
                    and np.isfinite(K_deep) and np.isfinite(K_shall)
                    and np.abs(K_deep - K_shall) < 1e6):
                dK_3pt = (K_deep - K_shall) / ds_total
            else:
                dK_3pt = np.nan

            # Eq. 8 terms
            t_vel    = 2.0 * K * dvds_si
            t_dK_sg  = vs_si * dK_sg
            t_dK_3pt = vs_si * dK_3pt

            pred_sg_pa  = eta * H * (t_vel + t_dK_sg)
            term_vel_pa = eta * H * t_vel

            # Eq. 11
            pred_eq11_pa = 0.1 * eta * K * vc_si

            records.append(dict(
                model        = model_label,
                meas_pa      = meas_pa,
                pred_sg_pa   = pred_sg_pa,
                pred_eq11_pa = pred_eq11_pa,
                term_vel_pa  = term_vel_pa,
                K            = K,
                H            = H,
                eta          = eta,
                vs_si        = vs_si,
                vc_si        = vc_si,
                Lv           = Lv,
                Lk_sg        = Lk_sg,
            ))

    return records


# ── Plotting helpers ──────────────────────────────────────────────────────────

K_THRESH = 0.0      # rad/m — exclude negative/zero curvature (unphysical slab states)

MPa = 1e6


def arr(records, key):
    return np.array([r[key] for r in records])


def model_colors(records):
    models = [r['model'] for r in records]
    unique = sorted(set(models))
    cmap   = matplotlib.colormaps.get_cmap('tab20').resampled(len(unique))
    color_of = {m: cmap(i) for i, m in enumerate(unique)}
    return models, unique, color_of


def corr_str(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if np.sum(ok) < 3:
        return "r=n/a"
    return f"N={np.sum(ok)}\nr={np.corrcoef(x[ok], y[ok])[0,1]:.2f}"


def one_to_one(ax, xdata, ydata, pad=0.05):
    both = np.concatenate([xdata[np.isfinite(xdata)], ydata[np.isfinite(ydata)]])
    if len(both) == 0:
        return
    lo, hi = np.nanpercentile(both, 1), np.nanpercentile(both, 99)
    span = hi - lo
    lo -= pad * span; hi += pad * span
    ax.plot([lo, hi], [lo, hi], 'k--', lw=0.9, zorder=0)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)


def label_box(ax, txt):
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va='top', fontsize=8,
            bbox=dict(fc='white', ec='none', alpha=0.75))


def _save(fig, out_dir, stem):
    for ext in ("png", "pdf"):
        path = os.path.join(out_dir, f"dqds_diag.{stem}.{ext}")
        fig.savefig(path, dpi=200, bbox_inches='tight')
        print(f"  saved: {path}")
    plt.close(fig)


# ── Figures ───────────────────────────────────────────────────────────────────

def make_figure_comparison(records, out_dir):
    """2-panel: Eq.11 uncalibrated then calibrated vs measured (per-model colored).
    Left shows systematic ~3x overprediction; right shows calibrated r=0.86 result."""
    models, unique, color_of = model_colors(records)

    meas_pa = arr(records, "meas_pa")
    eq11_pa = arr(records, "pred_eq11_pa")
    K_a     = arr(records, "K")

    # derive calibrated coefficient from K>0, same-sign subset
    ok_cal = (np.isfinite(meas_pa) & np.isfinite(eq11_pa)
              & (np.abs(eq11_pa) > 1e4)
              & (np.sign(meas_pa) == np.sign(eq11_pa))
              & (K_a > 0.0))
    med         = np.nanmedian(meas_pa[ok_cal] / eq11_pa[ok_cal])
    calib_coeff = 0.1 * med

    meas     = meas_pa    / MPa
    eq11     = eq11_pa    / MPa
    eq11_cal = eq11 * med

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    panels = [
        (axes[0], eq11,     "0.1 · η · K · v_c  [MPa]",
         "Eq. 11 (uncalibrated)"),
        (axes[1], eq11_cal, f"{calib_coeff:.4f} · η · K · v_c  [MPa]",
         f"Eq. 11 calibrated  (coeff = {calib_coeff:.4f})"),
    ]

    handles = []
    for ax, pred, xlabel, title in panels:
        for m in unique:
            idx = np.array([i for i, r in enumerate(records) if r['model'] == m])
            ok  = np.isfinite(pred[idx]) & np.isfinite(meas[idx])
            h = ax.scatter(pred[idx][ok], meas[idx][ok],
                           s=10, alpha=0.55, color=color_of[m],
                           edgecolors='none', label=m, zorder=3)
            if ax is axes[0]:
                handles.append(h)
        all_ok = np.isfinite(pred) & np.isfinite(meas)
        one_to_one(ax, pred[all_ok], meas[all_ok])
        ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Measured  dQ/ds  [MPa]")
        ax.set_title(title)
        ax.grid(True, ls='--', lw=0.4, alpha=0.4)
        label_box(ax, corr_str(pred, meas))

    fig.legend(handles, unique, title='Model', fontsize=6.5,
               loc='lower center', ncol=min(len(unique), 5),
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("dQ/ds scaling law  (paper sign convention)",
                 fontsize=10, fontweight='bold')
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    _save(fig, out_dir, "fig1_scaling_law")


def make_figure_calibration(records, out_dir):
    """3-panel: prefactor vs K (motivates stability filter),
    K-filtered prefactor histogram, calibrated Eq.11 vs measured."""
    models, unique, color_of = model_colors(records)

    meas  = arr(records, "meas_pa")
    eq11  = arr(records, "pred_eq11_pa")
    K_a   = arr(records, "K")

    # prefactor: same-sign, non-trivial eq11 only
    ok_base = (np.isfinite(meas) & np.isfinite(eq11)
               & (np.abs(eq11) > 1e4)
               & (np.sign(meas) == np.sign(eq11)))
    pf = np.where(ok_base, meas / eq11, np.nan)

    ok_filt = ok_base & (K_a > 0.0)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    c_pts = np.array([color_of[m] for m in models])

    # ── (a) prefactor vs K ───────────────────────────────────────────────
    ax = axes[0]
    valid = ok_base & np.isfinite(pf)
    ax.scatter(K_a[valid] * 1e6, pf[valid],
               s=10, alpha=0.5, c=c_pts[valid], edgecolors='none', zorder=3)
    ax.axvline(0.0, color='k', lw=1.4, ls='--', label='K = 0')
    ax.axhline(1.0, color='gray', lw=1.0, ls='--')
    ax.set_ylim(-0.2, np.nanpercentile(pf[valid], 97) * 1.3)
    ax.set_xlabel("K  [×10⁻⁶ rad/m]")
    ax.set_ylabel("meas  /  (0.1 η K v_c)")
    ax.set_title("Eq. 11 prefactor vs K\n(dashed = K > 0 filter)")
    ax.legend(fontsize=8); ax.grid(True, ls='--', lw=0.4, alpha=0.4)
    r_kpf = np.corrcoef(K_a[valid], pf[valid])[0, 1]
    label_box(ax, f"r = {r_kpf:.2f}")

    # ── (b) K-filtered prefactor histogram ──────────────────────────────
    ax = axes[1]
    v = pf[np.isfinite(pf) & ok_filt]
    p2, p98 = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    ax.hist(v[(v > p2) & (v < p98)], bins=35,
            color='steelblue', edgecolor='none', alpha=0.75)
    med = np.nanmedian(v)
    ax.axvline(med, color='k',    lw=1.8, label=f'median = {med:.3f}')
    ax.axvline(1.0, color='gray', lw=1.0, ls='--', label='1:1')
    ax.set_xlabel("meas  /  (0.1 η K v_c)")
    ax.set_ylabel("count")
    ax.set_title(f"Prefactor  (K > 0,  N = {len(v)})")
    ax.legend(fontsize=8); ax.grid(True, ls='--', lw=0.4, alpha=0.4)
    calib_coeff = 0.1 * med
    label_box(ax, f"median = {med:.3f}\n→ coeff = {calib_coeff:.4f}")

    # ── (c) calibrated scatter ───────────────────────────────────────────
    ax = axes[2]
    pred_calib = eq11 * med / MPa
    meas_mpa   = meas / MPa
    handles = []
    for m in unique:
        idx = np.array([i for i, r in enumerate(records) if r['model'] == m])
        ok  = ok_filt[idx] & np.isfinite(pred_calib[idx]) & np.isfinite(meas_mpa[idx])
        h = ax.scatter(pred_calib[idx][ok], meas_mpa[idx][ok],
                       s=12, alpha=0.6, color=color_of[m],
                       edgecolors='none', label=m, zorder=3)
        handles.append(h)
    all_ok = ok_filt & np.isfinite(pred_calib) & np.isfinite(meas_mpa)
    one_to_one(ax, pred_calib[all_ok], meas_mpa[all_ok])
    ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5)
    ax.set_xlabel(f"{calib_coeff:.4f} · η · K · v_c  [MPa]")
    ax.set_ylabel("Measured  dQ/ds  [MPa]")
    ax.set_title("Calibrated Eq. 11 vs measured\n(K-filtered)")
    ax.grid(True, ls='--', lw=0.4, alpha=0.5)
    label_box(ax, corr_str(pred_calib[ok_filt], meas_mpa[ok_filt]))

    fig.legend(handles=handles, title='Model', fontsize=6.5,
               loc='lower center', ncol=min(len(unique), 5),
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Eq. 11 calibration",
                 fontsize=10, fontweight='bold')
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    _save(fig, out_dir, "fig2_calibration")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    files = args if args else sorted(glob.glob("text_files/TESTC/*.txt"))
    if not files:
        print("No input files. Run from analysis/ or pass file paths.")
        sys.exit(1)

    out_dir = "plots/tmp"
    os.makedirs(out_dir, exist_ok=True)
    print(f"Loading {len(files)} file(s)...")
    records = load_records(files)
    print(f"  {len(records)} valid timestep records")
    if not records:
        print("No valid records found — check file paths and column count.")
        sys.exit(1)

    print("Making figures...")
    make_figure_comparison(records, out_dir)
    make_figure_calibration(records, out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
