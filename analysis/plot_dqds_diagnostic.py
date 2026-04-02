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
        dz_m = dz_km * 1e3

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

        for row in data:
            dip_rad = np.deg2rad(row[COL_DIP])
            if dip_rad <= 0 or not np.isfinite(dip_rad):
                continue

            H       = row[COL_H]
            K       = row[COL_K]
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
                model   = label,
                meas_pa = meas_pa,
                K       = K,
                H       = H,
                eta     = eta,
                vc_si   = vc_si,
                vs_si   = vs_si,
                L_eff   = L_eff,
                Lv      = row[COL_LV],
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


def corr_str(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return "r=n/a"
    return f"N={ok.sum()}\nr={np.corrcoef(x[ok], y[ok])[0,1]:.2f}"


def one_to_one(ax, x, y, pad=0.05):
    both = np.concatenate([x[np.isfinite(x)], y[np.isfinite(y)]])
    if len(both) == 0: return
    lo, hi = np.nanpercentile(both, 1), np.nanpercentile(both, 99)
    span = hi - lo
    lo -= pad * span; hi += pad * span
    ax.plot([lo, hi], [lo, hi], 'k--', lw=0.9, zorder=0)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)


def label_box(ax, txt):
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va='top', fontsize=8,
            bbox=dict(fc='white', ec='none', alpha=0.75))


def _save(fig, out_dir, stem):
    for ext in ('png', 'pdf'):
        p = os.path.join(out_dir, f"dqds_diag.{stem}.{ext}")
        fig.savefig(p, dpi=200, bbox_inches='tight')
        print(f"  saved: {p}")
    plt.close(fig)


# ── Figure 1: L_eff distribution ──────────────────────────────────────────────

def make_figure_Leff(records, out_dir):
    """3-panel: L_eff histogram, L_eff vs K (coloured by model), per-model box."""
    models, unique, color_of = model_colors(records)

    L_eff_km = arr(records, 'L_eff') / 1e3      # m → km
    K_um     = arr(records, 'K')     * 1e6       # rad/m → ×10⁻⁶
    c_pts    = np.array([color_of[m] for m in models])

    ok = np.isfinite(L_eff_km) & (L_eff_km > 0)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # ── (a) histogram ────────────────────────────────────────────────────
    ax = axes[0]
    v = L_eff_km[ok]
    p2, p98 = np.nanpercentile(v, 2), np.nanpercentile(v, 98)
    v_clip = v[(v > p2) & (v < p98) & (v < 10000)]
    ax.hist(v_clip, bins=40, color='steelblue', edgecolor='none', alpha=0.75)
    med = np.nanmedian(v)
    ax.axvline(med, color='k', lw=1.8, label=f'median = {med:.0f} km')
    ax.set_xlabel("L_eff  [km]")
    ax.set_ylabel("count")
    ax.set_xlim(0, 10000)
    ax.set_title(f"L_eff distribution  (K > 0,  N = {ok.sum()})")
    ax.legend(fontsize=8); ax.grid(True, ls='--', lw=0.4, alpha=0.4)
    label_box(ax, f"L_eff = ηHKv_c / dQ/ds\nmedian = {med:.0f} km")

    # ── (b) L_eff vs K with linear fit ───────────────────────────────────
    ax = axes[1]
    ax.scatter(K_um[ok], L_eff_km[ok],
               s=10, alpha=0.5, c=c_pts[ok], edgecolors='none', zorder=3)
    # linear fit through origin: L_eff = slope * K
    Kv  = (arr(records, 'K') * 1e6)[ok]          # ×10⁻⁶ rad/m
    Lv  = L_eff_km[ok]
    fin = np.isfinite(Kv) & np.isfinite(Lv) & (Lv < 10000)
    slope = np.sum(Kv[fin] * Lv[fin]) / np.sum(Kv[fin] ** 2)  # OLS through origin
    xfit  = np.linspace(0, K_um[ok].max(), 200)
    ax.plot(xfit, slope * xfit, 'k-', lw=1.5,
            label=f'fit: L_eff = {slope:.0f}·K  (×10⁶ m)')
    ax.axhline(med, color='k', lw=1.0, ls='--', alpha=0.4)
    ax.set_xlabel("K  [×10⁻⁶ rad/m]")
    ax.set_ylabel("L_eff  [km]")
    ax.set_title("L_eff vs K\n(if L_eff ∝ K → dQ/ds independent of K)")
    ax.set_ylim(0, 10000)
    r_val = np.corrcoef(Kv[fin], Lv[fin])[0, 1]
    label_box(ax, f"r = {r_val:.2f}\nslope = {slope:.0f} ×10⁶ m")
    ax.legend(fontsize=7); ax.grid(True, ls='--', lw=0.4, alpha=0.4)

    # ── (c) per-model L_eff — sorted by median, with IQR bars ───────────
    ax = axes[2]
    # collect per-model stats and sort by median
    model_stats = []
    for m in unique:
        idx  = np.array([j for j, r in enumerate(records) if r['model'] == m])
        vals = L_eff_km[idx]
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if len(vals) == 0:
            continue
        model_stats.append((np.median(vals), m, vals))
    model_stats.sort(key=lambda x: x[0])

    for i, (med_m, m, vals) in enumerate(model_stats):
        q25, q75 = np.percentile(vals, 25), np.percentile(vals, 75)
        col = color_of[m]
        # individual points (jittered)
        jitter = (np.random.default_rng(i).random(len(vals)) - 0.5) * 0.4
        ax.scatter(i + jitter, vals, s=8, alpha=0.35, color=col, edgecolors='none', zorder=2)
        # IQR bar
        ax.plot([i, i], [q25, q75], color=col, lw=3, solid_capstyle='round', zorder=3)
        # median tick
        ax.plot([i - 0.35, i + 0.35], [med_m, med_m],
                color=col, lw=2.5, solid_capstyle='round', zorder=4)

    ax.axhline(med, color='k', lw=1.2, ls='--', label=f'overall median = {med:.0f} km')
    sorted_labels = [m for _, m, _ in model_stats]
    ax.set_xticks(range(len(sorted_labels)))
    ax.set_xticklabels(sorted_labels, rotation=55, ha='right', fontsize=6)
    ax.set_ylabel("L_eff  [km]")
    ax.set_ylim(0, 10000)
    ax.set_title("L_eff per model  (sorted by median)\nbar = IQR, line = median")
    ax.legend(fontsize=8); ax.grid(True, ls='--', lw=0.4, alpha=0.4, axis='y')

    fig.suptitle("Effective length scale  L_eff = ηHKv_c / dQ/ds",
                 fontsize=10, fontweight='bold')
    fig.tight_layout()
    _save(fig, out_dir, "fig1_Leff_distribution")


# ── Figure 2: prediction using L_eff_median ───────────────────────────────────

def make_figure_prediction(records, out_dir):
    """2-panel: K-based prediction (all same-sign) vs no-K prediction (ηHvc/α).
    Tests whether K actually belongs in the scaling law."""
    models, unique, color_of = model_colors(records)

    meas_pa = arr(records, 'meas_pa')
    K_a     = arr(records, 'K')
    H_a     = arr(records, 'H')
    eta_a   = arr(records, 'eta')
    vc_a    = arr(records, 'vc_si')
    L_eff_a = arr(records, 'L_eff')

    # median L_eff (K>0, positive dQ/ds)
    ok_Leff = np.isfinite(L_eff_a) & (L_eff_a > 0)
    L_med   = np.nanmedian(L_eff_a[ok_Leff])

    # α = median(L_eff / K)  — slope of L_eff ∝ K fit (units m²)
    L_eff_km = L_eff_a / 1e3
    K_um     = K_a * 1e6
    fin      = ok_Leff & np.isfinite(K_um) & (L_eff_km < 10000)
    alpha    = np.sum(K_um[fin] * L_eff_km[fin]) / np.sum(K_um[fin] ** 2)
    # alpha in units (km / (×10⁻⁶ rad/m)) = km * m / (10⁻⁶) = 10⁶ km·m = 10⁹ m²
    # convert: α_m2 so that dQ/ds = η H vc / α_m2
    # L_eff [m] = alpha_scaled * K [m⁻¹]  → alpha_scaled [m²]
    alpha_m2 = (alpha * 1e3) / 1e-6    # km→m, ×10⁻⁶→m⁻¹

    pred_K_pa    = eta_a * H_a * K_a * vc_a / L_med    # with K
    pred_noK_pa  = eta_a * H_a * vc_a / alpha_m2        # no K

    meas_mpa     = meas_pa    / MPa
    pred_K_mpa   = pred_K_pa  / MPa
    pred_noK_mpa = pred_noK_pa / MPa

    ok_all_K   = (np.isfinite(meas_mpa) & np.isfinite(pred_K_mpa)
                  & (np.sign(meas_pa) == np.sign(pred_K_mpa)))
    ok_all_noK = (np.isfinite(meas_mpa) & np.isfinite(pred_noK_mpa)
                  & (np.sign(meas_pa) == np.sign(pred_noK_mpa)))

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    panels = [
        (axes[0], pred_K_mpa,   ok_all_K,
         f"η H K v_c / L_eff  [MPa]\n(L_eff = {L_med/1e3:.0f} km)"),
        (axes[1], pred_noK_mpa, ok_all_noK,
         f"η H v_c / α  [MPa]\n(α = {alpha_m2/1e12:.2f} ×10¹² m²,  no K)"),
    ]

    handles = []
    for ax, pred, mask, xlabel in panels:
        for m in unique:
            idx = np.array([i for i, r in enumerate(records) if r['model'] == m])
            ok  = mask[idx]
            h = ax.scatter(pred[idx][ok], meas_mpa[idx][ok],
                           s=10, alpha=0.55, color=color_of[m],
                           edgecolors='none', label=m, zorder=3)
            if ax is axes[0] and h not in handles:
                handles.append(h)
        one_to_one(ax, pred[mask], meas_mpa[mask])
        ax.axhline(0, color='0.7', lw=0.5); ax.axvline(0, color='0.7', lw=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Measured  dQ/ds  [MPa]")
        ax.grid(True, ls='--', lw=0.4, alpha=0.4)
        label_box(ax, corr_str(pred[mask], meas_mpa[mask]))

    fig.legend(handles, unique, title='Model', fontsize=6.5,
               loc='lower center', ncol=min(len(unique), 5),
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Does K belong in the dQ/ds scaling?  (panel 3: no-K formula)",
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

    L_vals = arr(records, 'L_eff')
    ok     = np.isfinite(L_vals) & (L_vals > 0)
    print(f"\n  L_eff  median = {np.nanmedian(L_vals[ok])/1e3:.0f} km"
          f"  (N={ok.sum()}, from K>0 and dQ/ds>0 timesteps)")

    print("\nMaking figures...")
    make_figure_Leff(records, out_dir)
    make_figure_prediction(records, out_dir)
    print("Done.")


if __name__ == "__main__":
    main()
