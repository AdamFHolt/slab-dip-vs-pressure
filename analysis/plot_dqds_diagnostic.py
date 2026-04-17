#!/usr/bin/env python3
"""
dQ/ds diagnostic script — L_eff framing.

The scaling law is expressed as:

    dQ/ds  =  η H K v_c / L_eff

where L_eff is the effective length scale implied by the data:

    L_eff  =  η H K v_c / |dQ/ds_meas|

Key questions:
  1. Is L_eff approximately constant across models/timesteps?
  2. How well does η H K v_c / L_eff_median predict measured dQ/ds?

L_eff is the quantity to carry forward to the Earth application section:
  dQ/ds / DP_anal  =  H K v_c / (L_eff · Δρ g cos θ)

Sign convention: paper's dQ/ds > 0 = force resisting slab-normal buoyancy.
  code col 6 (slab_stress_term) = -dQ/ds_paper, so meas_paper = -col6.

Overturned (rollover) slab models are tagged and shown as hollow symbols.
L_eff calibration uses non-overturned models only.

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
COL_DIP     = 5    # deg  — slab dip at analysis depth
COL_DQDS    = 6    # Pa   — slab_stress_term = -dQ/ds_paper
COL_H       = 9    # m    — slab-normal thickness
COL_K       = 11   # rad/m
COL_DK      = 12   # rad/m² — dK/ds
COL_K_SHALL = 13   # rad/m
COL_K_DEEP  = 14   # rad/m
COL_VC      = 19   # cm/yr
COL_VS      = 21   # m/yr
COL_DVDS    = 22   # yr⁻¹
COL_LV      = 23   # m
COL_LK      = 24   # m
COL_ETA     = 26   # Pa·s

SPY        = 365.25 * 24.0 * 3600.0
CMYR_TO_MS = 0.01 / SPY

# Models excluded from L_eff calibration — shown hollow
# Overturned (rollover) geometry:
OVERTURNED_MODELS = {
    'new_1000plates',         # η' = 1000, free plates
    'new_FixedOP_1000plates', # η' = 1000, fixed OP
    'FixedOP_lower-res_new',  # η' = 500,  fixed OP
    'new_FixedOP_375plates',  # η' = 375,  fixed OP
}

# Canonical model order: fixed SP / free / fixed OP, grouped by η'
MODEL_ORDER = [
    'new_FixedSP_50plates',    'new_50plates',    'new_FixedOP_50plates',
    'new_FixedSP_250plates',   'new_250plates',   'new_FixedOP_250plates',
    'new_FixedSP_375plates',   'new_375plates',   'new_FixedOP_375plates',
    'FixedSP_lower-res_new2',  'new2',            'FixedOP_lower-res_new',
    'new_FixedSP_1000plates2', 'new_1000plates',  'new_FixedOP_1000plates',
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def ensure_2d(a):
    return a.reshape(1, -1) if a.ndim == 1 else a


def parse_shear_dz_km(fpath):
    m = re.search(r'shear-dz([\d.]+?)\.(?:ds|km|[a-zA-Z])', os.path.basename(fpath))
    return float(m.group(1)) if m else None


# ── Data loading ──────────────────────────────────────────────────────────────

def load_records(files):
    records = []
    for fpath in files:
        dz_km = parse_shear_dz_km(fpath)
        if dz_km is None:
            print(f"  skip (no shear-dz token): {fpath}"); continue

        try:
            data = ensure_2d(np.loadtxt(fpath))
        except Exception as e:
            print(f"  skip (load error): {fpath} — {e}"); continue
        if data.shape[1] < 27:
            print(f"  skip ({data.shape[1]} cols): {fpath}"); continue

        label = os.path.basename(fpath).replace('.txt', '')
        for pfx in ('2D_compositional_subd_lower-res_', '2D_compositional_subd_'):
            if label.startswith(pfx):
                label = label[len(pfx):]; break
        label = re.sub(r'\.z[\d.]+\.shear.*', '', label)

        overturned = label in OVERTURNED_MODELS

        for row in data:
            dip_rad = np.deg2rad(row[COL_DIP])
            if dip_rad <= 0 or not np.isfinite(dip_rad):
                continue

            H       = row[COL_H]
            K       = row[COL_K]
            dKds    = row[COL_DK]
            eta     = row[COL_ETA]
            vc_si   = row[COL_VC]  * CMYR_TO_MS
            vs_si   = row[COL_VS]  / SPY
            meas_pa = -row[COL_DQDS]   # paper sign: +ve resists buoyancy

            # L_eff = η H K v_c / |dQ/ds|  (only where meas is positive and K > 0)
            denom = np.abs(meas_pa)
            if (K > 0 and np.isfinite(K) and np.isfinite(eta)
                    and np.isfinite(H) and np.isfinite(vc_si)
                    and denom > 1e3 and meas_pa > 0):
                L_eff = eta * H * K * vc_si / denom   # m
            else:
                L_eff = np.nan

            records.append(dict(
                model      = label,
                overturned = overturned,
                meas_pa    = meas_pa,
                K          = K,
                dKds       = dKds,
                H          = H,
                KH         = K * H,        # dimensionless thin-sheet parameter
                eta        = eta,
                vc_si      = vc_si,
                vs_si      = vs_si,
                L_eff      = L_eff,
                Lv         = row[COL_LV],
            ))

    return records


# ── Plotting helpers ───────────────────────────────────────────────────────────

MPa = 1e6


def arr(records, key):
    return np.array([r[key] for r in records])


def model_colors(records):
    models   = [r['model'] for r in records]
    unique   = sorted(set(models))
    cmap     = matplotlib.colormaps.get_cmap('tab20').resampled(len(unique))
    color_of = {m: cmap(i) for i, m in enumerate(unique)}
    return models, unique, color_of


def corr_str(x, y, label=''):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return "r=n/a"
    r = np.corrcoef(x[ok], y[ok])[0, 1]
    prefix = f"{label}\n" if label else ""
    return f"{prefix}N={ok.sum()},  r={r:.2f}"


def one_to_one(ax, x, y, pad=0.05):
    both = np.concatenate([x[np.isfinite(x)], y[np.isfinite(y)]])
    if len(both) == 0: return
    lo, hi = np.nanpercentile(both, 1), np.nanpercentile(both, 99)
    span = hi - lo
    lo -= pad * span; hi += pad * span
    ax.plot([lo, hi], [lo, hi], 'k--', lw=0.9, zorder=0)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)


def label_box(ax, txt):
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va='top', fontsize=7.5,
            bbox=dict(fc='white', ec='none', alpha=0.75))


def _save(fig, out_dir, stem):
    for ext in ('png', 'pdf'):
        p = os.path.join(out_dir, f"dqds_diag.{stem}.{ext}")
        fig.savefig(p, dpi=200, bbox_inches='tight')
        print(f"  saved: {p}")
    plt.close(fig)


def _scatter(ax, x, y, color, overturned, s=10, alpha=0.55):
    """Plot normal points filled, overturned points as hollow circles."""
    norm = ~overturned
    if norm.any():
        ax.scatter(x[norm], y[norm], s=s, alpha=alpha,
                   color=color[norm], edgecolors='none', zorder=3)
    if overturned.any():
        ax.scatter(x[overturned], y[overturned], s=s+4, alpha=alpha+0.1,
                   facecolors='none', edgecolors=color[overturned],
                   linewidths=0.8, zorder=4)


# ── Figure 1: L_eff distribution ──────────────────────────────────────────────

def make_figure_Leff(records, out_dir):
    """3-panel: L_eff histogram, L_eff vs K, per-model box.
    Overturned models shown as hollow symbols; calibration uses normal models only."""
    models, unique, color_of = model_colors(records)

    L_eff_km  = arr(records, 'L_eff') / 1e3
    K_um      = arr(records, 'K')     * 1e6
    is_over   = arr(records, 'overturned').astype(bool)
    c_pts     = np.array([color_of[m] for m in models])

    ok      = np.isfinite(L_eff_km) & (L_eff_km > 0)
    ok_norm = ok & ~is_over
    ok_over = ok &  is_over

    med_all  = np.nanmedian(L_eff_km[ok])
    med_norm = np.nanmedian(L_eff_km[ok_norm])

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # ── (a) histogram ────────────────────────────────────────────────────
    ax = axes[0]
    clip = lambda v: v[(v > np.nanpercentile(v, 2)) & (v < 10000)]
    ax.hist(clip(L_eff_km[ok_norm]), bins=35, color='steelblue',
            edgecolor='none', alpha=0.75, label=f'normal  (median = {med_norm:.0f} km)')
    ax.hist(clip(L_eff_km[ok_over]), bins=15, histtype='step',
            edgecolor='tomato', linewidth=1.5, alpha=0.9,
            label=f'overturned  (median = {med_all:.0f} km all)')
    ax.axvline(med_norm, color='steelblue', lw=1.8, ls='-')
    ax.axvline(med_all,  color='k',         lw=1.2, ls='--')
    ax.set_xlabel("L_eff  [km]")
    ax.set_ylabel("count")
    ax.set_xlim(0, 10000)
    ax.set_title(f"L_eff distribution\nfilled: normal (N={ok_norm.sum()}), "
                 f"hollow: overturned (N={ok_over.sum()})")
    ax.legend(fontsize=7); ax.grid(True, ls='--', lw=0.4, alpha=0.4)

    # ── (b) L_eff vs K with linear fit (normal only) ─────────────────────
    ax = axes[1]
    _scatter(ax, K_um[ok], L_eff_km[ok], c_pts[ok], is_over[ok])
    # OLS fit through origin on normal models only
    Kv  = K_um[ok_norm]
    Lv  = L_eff_km[ok_norm]
    fin = np.isfinite(Kv) & np.isfinite(Lv) & (Lv < 10000)
    slope = np.sum(Kv[fin] * Lv[fin]) / np.sum(Kv[fin] ** 2)
    xfit  = np.linspace(0, K_um[ok].max(), 200)
    ax.plot(xfit, slope * xfit, 'k-', lw=1.5,
            label=f'fit (normal): L_eff = {slope:.0f}·K  (×10⁶ m)')
    r_val = np.corrcoef(Kv[fin], Lv[fin])[0, 1]
    label_box(ax, f"r = {r_val:.2f}  (normal only)\nslope = {slope:.0f} ×10⁶ m")
    ax.set_xlabel("K  [×10⁻⁶ rad/m]")
    ax.set_ylabel("L_eff  [km]")
    ax.set_title("L_eff vs K\n(hollow = overturned)")
    ax.set_ylim(0, 10000)
    ax.legend(fontsize=7); ax.grid(True, ls='--', lw=0.4, alpha=0.4)

    # ── (c) per-model L_eff — sorted by median, with IQR bars ───────────
    ax = axes[2]
    # order by canonical model sequence (η' group × BC type)
    present = set(unique)
    ordered = [m for m in MODEL_ORDER if m in present]
    ordered += sorted(present - set(MODEL_ORDER))   # any unlisted models appended

    model_stats = []
    for m in ordered:
        idx  = np.array([j for j, r in enumerate(records) if r['model'] == m])
        vals = L_eff_km[idx]
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if len(vals) == 0:
            continue
        is_ov = m in OVERTURNED_MODELS
        model_stats.append((np.median(vals), m, vals, is_ov))

    for i, (med_m, m, vals, is_ov) in enumerate(model_stats):
        q25, q75 = np.percentile(vals, 25), np.percentile(vals, 75)
        col = color_of[m]
        jitter = (np.random.default_rng(i).random(len(vals)) - 0.5) * 0.4
        if is_ov:
            ax.scatter(i + jitter, vals, s=8, alpha=0.35,
                       facecolors='none', edgecolors=col, linewidths=0.6, zorder=2)
            ax.plot([i, i], [q25, q75], color=col, lw=2, ls='--',
                    solid_capstyle='round', zorder=3)
            ax.plot([i - 0.35, i + 0.35], [med_m, med_m],
                    color=col, lw=2, ls='--', solid_capstyle='round', zorder=4)
        else:
            ax.scatter(i + jitter, vals, s=8, alpha=0.35,
                       color=col, edgecolors='none', zorder=2)
            ax.plot([i, i], [q25, q75], color=col, lw=3,
                    solid_capstyle='round', zorder=3)
            ax.plot([i - 0.35, i + 0.35], [med_m, med_m],
                    color=col, lw=2.5, solid_capstyle='round', zorder=4)

    ax.axhline(med_norm, color='steelblue', lw=1.2, ls='-',
               label=f'normal median = {med_norm:.0f} km')
    ax.axhline(med_all,  color='k',         lw=1.2, ls='--',
               label=f'all median = {med_all:.0f} km')
    sorted_labels = [m for _, m, _, _ in model_stats]
    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(sorted_labels, rotation=55, ha='right', fontsize=6)
    ax.set_ylabel("L_eff  [km]")
    ax.set_ylim(0, 10000)
    ax.set_title("L_eff per model  (fixedSP / free / fixedOP  ×  η')\ndashed = overturned")
    ax.legend(fontsize=7); ax.grid(True, ls='--', lw=0.4, alpha=0.4, axis='y')

    fig.suptitle("Effective length scale  L_eff = ηHKv_c / dQ/ds",
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    _save(fig, out_dir, "fig1_Leff_distribution")


# ── Figure 2: prediction using L_eff_median ───────────────────────────────────

def make_figure_prediction(records, out_dir):
    """3-panel: K-based (L_eff), no-K (α), and direct ηHvs·dK/ds prediction.
    L_eff calibrated from normal (non-overturned) models only.
    Overturned points shown as hollow circles."""
    models, unique, color_of = model_colors(records)

    meas_pa  = arr(records, 'meas_pa')
    K_a      = arr(records, 'K')
    dKds_a   = arr(records, 'dKds')
    H_a      = arr(records, 'H')
    eta_a    = arr(records, 'eta')
    vc_a     = arr(records, 'vc_si')
    vs_a     = arr(records, 'vs_si')
    L_eff_a  = arr(records, 'L_eff')
    is_over  = arr(records, 'overturned').astype(bool)
    c_pts    = np.array([color_of[m] for m in models])

    # calibrate L_eff and α from normal models only
    ok_Leff      = np.isfinite(L_eff_a) & (L_eff_a > 0)
    ok_Leff_norm = ok_Leff & ~is_over
    L_med_norm   = np.nanmedian(L_eff_a[ok_Leff_norm])
    L_med_all    = np.nanmedian(L_eff_a[ok_Leff])

    L_eff_km = L_eff_a / 1e3
    K_um     = K_a * 1e6
    fin      = ok_Leff_norm & np.isfinite(K_um) & (L_eff_km < 10000)
    alpha    = np.sum(K_um[fin] * L_eff_km[fin]) / np.sum(K_um[fin] ** 2)
    alpha_m2 = (alpha * 1e3) / 1e-6

    pred_K_pa    = eta_a * H_a * K_a * vc_a / L_med_norm
    pred_noK_pa  = eta_a * H_a * vc_a / alpha_m2
    pred_dKds_pa = eta_a * H_a * vs_a * dKds_a

    meas_mpa      = meas_pa      / MPa
    pred_K_mpa    = pred_K_pa    / MPa
    pred_noK_mpa  = pred_noK_pa  / MPa
    pred_dKds_mpa = pred_dKds_pa / MPa

    KH_a    = arr(records, 'KH')
    ok_K    = np.isfinite(meas_mpa) & np.isfinite(pred_K_mpa)    & (np.sign(meas_pa) == np.sign(pred_K_pa))
    ok_noK  = np.isfinite(meas_mpa) & np.isfinite(pred_noK_mpa)  & (np.sign(meas_pa) == np.sign(pred_noK_pa))
    ok_dK   = np.isfinite(meas_mpa) & np.isfinite(pred_dKds_mpa) & (np.sign(meas_pa) == np.sign(pred_dKds_pa))
    ok_dK_thin = ok_dK & (KH_a < 0.1)    # strict thin-sheet filter
    ok_dK_mod  = ok_dK & (KH_a < 0.2)    # moderate thin-sheet filter

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    panels = [
        (axes[0], pred_K_mpa,    ok_K,
         f"η H K v_c / L_eff  [MPa]\n(L_eff = {L_med_norm/1e3:.0f} km, normal only)"),
        (axes[1], pred_noK_mpa,  ok_noK,
         f"η H v_c / α  [MPa]\n(α = {alpha_m2/1e12:.2f} ×10¹² m²,  no K)"),
        (axes[2], pred_dKds_mpa, ok_dK,
         "η H v_s · dK/ds  [MPa]\n(dominant term, Eq. 8)"),
    ]

    legend_handles = []
    for ax, pred, mask, xlabel in panels:
        for m in unique:
            idx  = np.array([i for i, r in enumerate(records) if r['model'] == m])
            ok_m = mask[idx]
            col  = np.array([color_of[m]] * ok_m.sum())
            ov   = is_over[idx][ok_m]
            x    = pred[idx][ok_m]
            y    = meas_mpa[idx][ok_m]
            _scatter(ax, x, y, col, ov)
            if ax is axes[0]:
                # one proxy handle per model for legend
                h = ax.scatter([], [], s=10, color=color_of[m], label=m)
                legend_handles.append(h)

        if ax is axes[2]:
            m_lo = np.nanpercentile(meas_mpa[mask], 1)
            m_hi = np.nanpercentile(meas_mpa[mask], 99)
            span = m_hi - m_lo
            ax.set_xlim(m_lo - 0.05*span, m_hi + 0.05*span)
            ax.set_ylim(m_lo - 0.05*span, m_hi + 0.05*span)
            ax.plot([m_lo, m_hi], [m_lo, m_hi], 'k--', lw=0.9, zorder=0)
        else:
            one_to_one(ax, pred[mask], meas_mpa[mask])

        ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Measured  dQ/ds  [MPa]")
        ax.grid(True, ls='--', lw=0.4, alpha=0.4)

        # stats for normal and all subsets; extra KH filters for dK/ds panel
        norm_mask = mask & ~is_over
        all_mask  = mask
        s_norm = corr_str(pred[norm_mask], meas_mpa[norm_mask], 'excl. overturned')
        s_all  = corr_str(pred[all_mask],  meas_mpa[all_mask],  'all')
        if ax is axes[2]:
            s_thin = corr_str(pred[ok_dK_thin], meas_mpa[ok_dK_thin], 'KH < 0.1')
            s_mod  = corr_str(pred[ok_dK_mod],  meas_mpa[ok_dK_mod],  'KH < 0.2')
            label_box(ax, s_all + '\n' + s_mod + '\n' + s_thin)
        else:
            label_box(ax, s_norm + '\n' + s_all)

    # hollow circle proxy for overturned
    h_over = ax.scatter([], [], s=12, facecolors='none', edgecolors='0.3',
                        linewidths=0.8, label='overturned')
    fig.legend(legend_handles + [h_over],
               [m for m in unique] + ['overturned'],
               title='Model', fontsize=6.5,
               loc='lower center', ncol=min(len(unique) + 1, 6),
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"dQ/ds predictions  (L_eff calibrated on normal slabs: "
                 f"{L_med_norm/1e3:.0f} km;  all: {L_med_all/1e3:.0f} km)",
                 fontsize=10, fontweight='bold')
    fig.tight_layout(rect=[0, 0.10, 1, 1])
    _save(fig, out_dir, "fig2_Leff_prediction")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args  = [a for a in sys.argv[1:] if not a.startswith('--')]
    files = args if args else sorted(glob.glob('text_files/TESTC/*.txt'))
    if not files:
        print("No input files. Run from analysis/ or pass file paths.")
        sys.exit(1)

    out_dir = 'plots/tmp'
    os.makedirs(out_dir, exist_ok=True)
    print(f"Loading {len(files)} file(s)...")
    records = load_records(files)
    print(f"  {len(records)} valid records")
    if not records:
        print("No valid records — check paths and column count.")
        sys.exit(1)

    is_over = arr(records, 'overturned').astype(bool)
    for label, mask in [('all models', np.ones(len(records), bool)),
                        ('excl. overturned', ~is_over)]:
        L_vals = arr(records, 'L_eff')[mask]
        ok = np.isfinite(L_vals) & (L_vals > 0)
        print(f"  L_eff median = {np.nanmedian(L_vals[ok])/1e3:.0f} km"
              f"  (N={ok.sum()})  [{label}]")

    print("\nMaking figures...")
    make_figure_Leff(records, out_dir)
    make_figure_prediction(records, out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
