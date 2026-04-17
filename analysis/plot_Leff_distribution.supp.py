#!/bin/python
"""
Supplementary figure: L_eff distribution (panels 1 and 3 from diagnostic fig1).

  Panel 1: L_eff histogram — normal (filled) vs overturned (outlined)
  Panel 2: L_eff per model — median marks, canonical model order

Auto-loads text_files/TESTC/*.txt.
Output: plots/DP-comparisons/compilations/Leff_distribution.supp.{png,pdf}

Usage:
    python3 plot_Leff_distribution.supp.py
"""

import os
import re
import glob
import numpy as np
import matplotlib
import matplotlib as mpl
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.font_manager as fm

font_path = "/home/holt/.local/share/fonts/MYRIADPRO-REGULAR.OTF"
myriad_pro = fm.FontProperties(fname=font_path)

mpl.rcParams['font.family'] = 'Myriad Pro'
mpl.rcParams['font.size'] = 9
mpl.rcParams['axes.labelsize'] = 9
mpl.rcParams['axes.labelpad'] = 1.5
mpl.rcParams['xtick.labelsize'] = 7.5
mpl.rcParams['ytick.labelsize'] = 7.5
mpl.rcParams['xtick.major.pad'] = 2
mpl.rcParams['ytick.major.pad'] = 2
mpl.rcParams['xtick.major.size'] = 2.5
mpl.rcParams['ytick.major.size'] = 2.5
mpl.rcParams['xtick.minor.size'] = 1.25
mpl.rcParams['ytick.minor.size'] = 1.25

# ── Column indices ─────────────────────────────────────────────────────────────
COL_DIP  = 5
COL_DQDS = 6
COL_H    = 9
COL_K    = 11
COL_VC   = 19
COL_ETA  = 26

SPY        = 365.25 * 24.0 * 3600.0
CMYR_TO_MS = 0.01 / SPY

OVERTURNED_MODELS = {
    'new_1000plates',
    'new_FixedOP_1000plates',
    'FixedOP_lower-res_new',
    'new_FixedOP_375plates',  # η' = 375, fixed OP
}

MODEL_ORDER = [
    'new_FixedSP_50plates',    'new_50plates',    'new_FixedOP_50plates',
    'new_FixedSP_250plates',   'new_250plates',   'new_FixedOP_250plates',
    'new_FixedSP_375plates',   'new_375plates',   'new_FixedOP_375plates',
    'FixedSP_lower-res_new2',  'new2',            'FixedOP_lower-res_new',
    'new_FixedSP_1000plates2', 'new_1000plates',  'new_FixedOP_1000plates',
]

MODEL_LABELS = {
    'new_FixedSP_50plates':    r"fixed SP, $\eta'$=50",
    'new_50plates':            r"free, $\eta'$=50",
    'new_FixedOP_50plates':    r"fixed OP, $\eta'$=50",
    'new_FixedSP_250plates':   r"fixed SP, $\eta'$=250",
    'new_250plates':           r"free, $\eta'$=250",
    'new_FixedOP_250plates':   r"fixed OP, $\eta'$=250",
    'new_FixedSP_375plates':   r"fixed SP, $\eta'$=375",
    'new_375plates':           r"free, $\eta'$=375",
    'new_FixedOP_375plates':   r"fixed OP, $\eta'$=375",
    'FixedSP_lower-res_new2':  r"fixed SP, $\eta'$=500",
    'new2':                    r"free, $\eta'$=500",
    'FixedOP_lower-res_new':   r"fixed OP, $\eta'$=500",
    'new_FixedSP_1000plates2': r"fixed SP, $\eta'$=1000",
    'new_1000plates':          r"free, $\eta'$=1000",
    'new_FixedOP_1000plates':  r"fixed OP, $\eta'$=1000",
}

# ── Data loading ───────────────────────────────────────────────────────────────

def ensure_2d(a):
    return a.reshape(1, -1) if a.ndim == 1 else a


def load_records(files):
    records = []
    for fpath in files:
        try:
            data = ensure_2d(np.loadtxt(fpath))
        except Exception:
            continue
        if data.shape[1] < 27:
            continue

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
            H, K, eta = row[COL_H], row[COL_K], row[COL_ETA]
            vc_si   = row[COL_VC] * CMYR_TO_MS
            meas_pa = -row[COL_DQDS]
            denom   = np.abs(meas_pa)
            if (K > 0 and np.isfinite(K) and np.isfinite(eta)
                    and np.isfinite(H) and np.isfinite(vc_si)
                    and denom > 1e3 and meas_pa > 0):
                L_eff = eta * H * K * vc_si / denom
            else:
                L_eff = np.nan
            records.append(dict(model=label, overturned=overturned, L_eff=L_eff))
    return records


# ── Plot ───────────────────────────────────────────────────────────────────────

def make_figure(records):
    models   = [r['model'] for r in records]
    unique   = sorted(set(models))
    cmap     = matplotlib.colormaps.get_cmap('tab20').resampled(len(unique))
    color_of = {m: cmap(i) for i, m in enumerate(unique)}

    L_eff_km = np.array([r['L_eff'] for r in records]) / 1e3
    is_over  = np.array([r['overturned'] for r in records], dtype=bool)

    ok      = np.isfinite(L_eff_km) & (L_eff_km > 0)
    ok_norm = ok & ~is_over
    ok_over = ok &  is_over

    med_norm = np.nanmedian(L_eff_km[ok_norm])
    med_all  = np.nanmedian(L_eff_km[ok])

    fig = plt.figure(figsize=(10, 4))
    gs  = GridSpec(1, 2, figure=fig, wspace=0.38)

    # ── Panel 1: histogram ────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])

    clip = lambda v: v[(v > np.nanpercentile(v, 2)) & (v < 10000)]
    ax1.hist(clip(L_eff_km[ok]), bins=35, color='gray',
             edgecolor='none', alpha=0.5,
             label=f'all  (N={ok.sum()}, median={med_all:.0f} km)')
    ax1.hist(clip(L_eff_km[ok_norm]), bins=35, color='steelblue',
             edgecolor='none', alpha=0.75,
             label=f'non-rollover  (N={ok_norm.sum()}, median={med_norm:.0f} km)')
    ax1.axvline(med_all,  color='k',         lw=1.0, ls='--')
    ax1.axvline(med_norm, color='steelblue', lw=1.5, ls='-')
    ax1.set_xlabel(r'$L_\mathrm{eff}$  [km]')
    ax1.set_ylabel('count')
    ax1.set_xlim(0, 7000)
    ax1.xaxis.set_minor_locator(plt.MultipleLocator(500))
    ax1.yaxis.set_minor_locator(plt.MultipleLocator(5))
    ax1.grid(True, which='major', color='lightgray', linestyle='--', linewidth=0.5, zorder=0)
    ax1.legend(fontsize=7.5, framealpha=0.8)

    # ── Panel 2: per-model strip ───────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])

    present = set(unique)
    ordered = [m for m in MODEL_ORDER if m in present]
    ordered += sorted(present - set(MODEL_ORDER))

    model_stats = []
    for m in ordered:
        vals = L_eff_km[[j for j, r in enumerate(records) if r['model'] == m]]
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if len(vals) == 0:
            continue
        is_ov = m in OVERTURNED_MODELS
        model_stats.append((np.median(vals), m, vals, is_ov))

    for i, (med_m, m, vals, is_ov) in enumerate(model_stats):
        col = color_of[m]
        rng = np.random.default_rng(i)
        jitter = (rng.random(len(vals)) - 0.5) * 0.4
        if is_ov:
            ax2.scatter(i + jitter, vals, s=5, alpha=0.3,
                        facecolors='none', edgecolors=col, linewidths=0.5, zorder=2)
            ax2.plot([i - 0.35, i + 0.35], [med_m, med_m],
                     color=col, lw=1.5, ls='--', zorder=4)
        else:
            ax2.scatter(i + jitter, vals, s=5, alpha=0.3,
                        color=col, edgecolors='none', zorder=2)
            ax2.plot([i - 0.35, i + 0.35], [med_m, med_m],
                     color=col, lw=2.0, zorder=4)

    ax2.axhline(med_all,  color='k',         lw=1.0, ls='--',
                label=f'all median = {med_all:.0f} km')
    ax2.axhline(med_norm, color='steelblue', lw=1.5, ls='-',
                label=f'non-rollover median = {med_norm:.0f} km')

    tick_labels = [MODEL_LABELS.get(m, m) for _, m, _, _ in model_stats]
    ax2.set_xticks(range(len(tick_labels)))
    ax2.set_xticklabels(tick_labels, rotation=55, ha='right', fontsize=7)
    ax2.set_ylabel(r'$L_\mathrm{eff}$  [km]')
    ax2.set_ylim(0, 7000)
    ax2.yaxis.set_minor_locator(plt.MultipleLocator(500))
    ax2.grid(True, which='major', color='lightgray', linestyle='--',
             linewidth=0.5, zorder=0, axis='y')
    ax2.legend(fontsize=7.5, framealpha=0.8)

    return fig


# ── Main ───────────────────────────────────────────────────────────────────────

files = sorted(glob.glob('text_files/TESTC/*.txt'))
if not files:
    print("No files found in text_files/TESTC/. Run from analysis/.")
    raise SystemExit(1)

print(f"Loading {len(files)} file(s)...")
records = load_records(files)
print(f"  {len(records)} valid records")

out_dir = 'plots/DP-comparisons/compilations'
os.makedirs(out_dir, exist_ok=True)

fig = make_figure(records)
for ext in ('png', 'pdf'):
    p = os.path.join(out_dir, f'Leff_distribution.supp.{ext}')
    fig.savefig(p, dpi=600, bbox_inches='tight', format=ext)
    print(f"  saved: {p}")
plt.close(fig)
print("Done.")
