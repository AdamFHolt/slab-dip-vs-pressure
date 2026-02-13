#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Columns in extract_properties.py output
COL_TIME = 0
COL_DQDS_MEASURED = 6
COL_DQDS_SPLITL = 25


def ensure_2d(arr):
    if arr.ndim == 1:
        return arr.reshape(1, -1)
    return arr


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_dqds_splitL.py <text_file_1> [<text_file_2> ...]")
        sys.exit(1)

    files = sys.argv[1:]
    out_dir = "plots/DP-comparisons/compilations"
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax_ts, ax_sc = axes

    all_meas = []
    all_pred = []

    for fpath in files:
        data = np.loadtxt(fpath)
        data = ensure_2d(data)
        base = os.path.basename(fpath).replace(".txt", "")

        if data.shape[1] <= COL_DQDS_SPLITL:
            print(f"Skipping {fpath}: missing split-L columns (need >= {COL_DQDS_SPLITL + 1} columns).")
            continue

        t = data[:, COL_TIME]
        meas = data[:, COL_DQDS_MEASURED] / 1.0e6
        pred = data[:, COL_DQDS_SPLITL] / 1.0e6

        finite = np.isfinite(meas) & np.isfinite(pred)
        t = t[finite]
        meas = meas[finite]
        pred = pred[finite]

        if len(t) == 0:
            print(f"Skipping {fpath}: no finite measured/predicted dQ/ds values.")
            continue

        all_meas.append(meas)
        all_pred.append(pred)

        ax_ts.plot(t, meas, "-", linewidth=1.5, label=f"{base} measured")
        ax_ts.plot(t, pred, "--", linewidth=1.5, label=f"{base} split-L")

        ax_sc.scatter(meas, pred, s=16, alpha=0.7, label=base)

    if not all_meas:
        print("No valid files to plot.")
        sys.exit(2)

    meas_cat = np.concatenate(all_meas)
    pred_cat = np.concatenate(all_pred)

    lim_min = min(np.min(meas_cat), np.min(pred_cat))
    lim_max = max(np.max(meas_cat), np.max(pred_cat))
    if lim_max > lim_min:
        ax_sc.plot([lim_min, lim_max], [lim_min, lim_max], "k--", linewidth=1.0)

    if len(meas_cat) > 1:
        corr = np.corrcoef(meas_cat, pred_cat)[0, 1]
        ax_sc.text(
            0.03,
            0.97,
            f"N={len(meas_cat)}\nr={corr:.2f}",
            transform=ax_sc.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8),
        )

    ax_ts.set_xlabel("Timestep")
    ax_ts.set_ylabel("dQ/ds [MPa]")
    ax_ts.set_title("Measured vs Split-L dQ/ds (time)")
    ax_ts.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax_ts.legend(fontsize=7, loc="best")

    ax_sc.set_xlabel("Measured dQ/ds [MPa]")
    ax_sc.set_ylabel("Predicted dQ/ds [MPa]")
    ax_sc.set_title("Split-L vs Measured dQ/ds")
    ax_sc.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax_sc.legend(fontsize=7, loc="best")

    plt.tight_layout()

    if len(files) == 1:
        stem = os.path.basename(files[0]).replace(".txt", "")
    else:
        stem = "multi-model"
    out_png = os.path.join(out_dir, f"dqds_splitL_compare.{stem}.png")
    out_pdf = os.path.join(out_dir, f"dqds_splitL_compare.{stem}.pdf")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
