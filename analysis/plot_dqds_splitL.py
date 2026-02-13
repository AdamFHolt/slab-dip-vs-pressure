#!/usr/bin/env python3
import os
import sys
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Columns in extract_properties.py output
COL_TIME = 0
COL_DQDS_MEASURED = 6
COL_H = 9
COL_K = 11
COL_LV = 23
COL_LK = 24
COL_DQDS_SPLITL = 25
COL_ETA = 26
COL_VS = 21
COL_L_TOTAL_INV = 29
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0


def ensure_2d(arr):
    if arr.ndim == 1:
        return arr.reshape(1, -1)
    return arr


def main():
    args = sys.argv[1:]
    use_lv_as_full = False
    use_2lv_as_full = False
    const_l_km = None
    if "--full-L-equals-Lv" in args:
        use_lv_as_full = True
        args = [a for a in args if a != "--full-L-equals-Lv"]
    if "--full-L-equals-2Lv" in args:
        use_2lv_as_full = True
        args = [a for a in args if a != "--full-L-equals-2Lv"]
    if "--const-L-km" in args:
        idx = args.index("--const-L-km")
        if idx + 1 >= len(args):
            print("Missing value for --const-L-km")
            sys.exit(1)
        try:
            const_l_km = float(args[idx + 1])
        except ValueError:
            print("Invalid value for --const-L-km; expected a number in km.")
            sys.exit(1)
        if const_l_km <= 0:
            print("--const-L-km must be > 0")
            sys.exit(1)
        args = args[:idx] + args[idx + 2 :]
    mode_count = int(use_lv_as_full) + int(use_2lv_as_full) + int(const_l_km is not None)
    if mode_count > 1:
        print("Use only one mode: --full-L-equals-Lv, --full-L-equals-2Lv, or --const-L-km <value>")
        sys.exit(1)

    if len(args) < 1:
        files = sorted(glob.glob("text_files/TESTC/*.txt"))
        if not files:
            print("No input files provided and none found in text_files/TESTC/*.txt")
            print("Usage: python3 plot_dqds_splitL.py [--full-L-equals-Lv] [--full-L-equals-2Lv] [--const-L-km <value>] <text_file_1> [<text_file_2> ...]")
            sys.exit(1)
        print(f"Using all TESTC files: {len(files)}")
    else:
        files = args
    out_dir = "plots/DP-comparisons/compilations"
    os.makedirs(out_dir, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    ax_sc = axes[0, 0]
    ax_lv = axes[0, 1]
    ax_lk = axes[1, 0]
    ax_lt = axes[1, 1]

    all_meas = []
    all_pred = []

    for fpath in files:
        data = np.loadtxt(fpath)
        data = ensure_2d(data)
        base = os.path.basename(fpath).replace(".txt", "")

        if data.shape[1] <= COL_L_TOTAL_INV:
            print(f"Skipping {fpath}: missing columns (need >= {COL_L_TOTAL_INV + 1} columns).")
            continue

        t = data[:, COL_TIME]
        meas = -data[:, COL_DQDS_MEASURED] / 1.0e6
        h = data[:, COL_H]
        k = data[:, COL_K]
        vs = data[:, COL_VS]
        lv = data[:, COL_LV]
        lk = data[:, COL_LK]
        eta = data[:, COL_ETA]
        lt_inv_saved = data[:, COL_L_TOTAL_INV]

        h_over_lv = np.full_like(lv, np.nan, dtype=float)
        h_over_lk = np.full_like(lk, np.nan, dtype=float)
        h_over_l = np.full_like(lt_inv_saved, np.nan, dtype=float)
        valid_h = np.isfinite(h) & (h > 0)
        valid_lv = np.isfinite(lv) & (lv > 0)
        valid_lk = np.isfinite(lk) & (lk > 0)
        valid_lt_saved = np.isfinite(lt_inv_saved) & (lt_inv_saved > 0)
        lt_inv = np.full_like(lt_inv_saved, np.nan, dtype=float)
        if const_l_km is not None:
            lt_inv[:] = 1.0 / (const_l_km * 1.0e3)
        elif use_2lv_as_full:
            lt_inv[valid_lv] = 2.0 / lv[valid_lv]
        elif use_lv_as_full:
            lt_inv[valid_lv] = 1.0 / lv[valid_lv]
        else:
            lt_inv[valid_lt_saved] = lt_inv_saved[valid_lt_saved]

        # Recompute predicted dQ/ds using chosen full-length term.
        # vs from model output is m/yr, so convert to m/s for SI-consistent stress.
        pred_pa = eta * h * (k * (vs / SECONDS_PER_YEAR)) * lt_inv
        pred = pred_pa / 1.0e6

        valid_lt = np.isfinite(lt_inv) & (lt_inv > 0)
        valid_h_lv = valid_h & valid_lv
        valid_h_lk = valid_h & valid_lk
        valid_h_lt = valid_h & valid_lt
        h_over_lv[valid_h_lv] = h[valid_h_lv] / lv[valid_h_lv]
        h_over_lk[valid_h_lk] = h[valid_h_lk] / lk[valid_h_lk]
        h_over_l[valid_h_lt] = h[valid_h_lt] * lt_inv[valid_h_lt]

        finite_pred = np.isfinite(meas) & np.isfinite(pred)
        if np.any(finite_pred):
            all_meas.append(meas[finite_pred])
            all_pred.append(pred[finite_pred])
            ax_sc.scatter(pred[finite_pred], meas[finite_pred], s=16, alpha=0.7, label=base)
        else:
            print(f"Skipping dQ/ds scatter for {fpath}: no finite measured/predicted values.")

        finite_lv = np.isfinite(t) & np.isfinite(h_over_lv)
        if np.any(finite_lv):
            ax_lv.plot(t[finite_lv], h_over_lv[finite_lv], linewidth=1.3, label=base)

        finite_lk = np.isfinite(t) & np.isfinite(h_over_lk)
        if np.any(finite_lk):
            ax_lk.plot(t[finite_lk], h_over_lk[finite_lk], linewidth=1.3, label=base)

        finite_hlt = np.isfinite(t) & np.isfinite(h_over_l)
        if np.any(finite_hlt):
            ax_lt.plot(t[finite_hlt], h_over_l[finite_hlt], linewidth=1.3, label=base)

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

    ax_sc.set_xlabel("Predicted dQ/ds [MPa]")
    ax_sc.set_ylabel("-Measured dQ/ds [MPa]")
    if const_l_km is not None:
        ax_sc.set_title(f"Pred vs Measured dQ/ds (L_full = {const_l_km:g} km)")
    elif use_2lv_as_full:
        ax_sc.set_title("Pred vs Measured dQ/ds (L_full = L_v/2; 1/L = 2/L_v)")
    elif use_lv_as_full:
        ax_sc.set_title("Pred vs Measured dQ/ds (L_full = L_v)")
    else:
        ax_sc.set_title("Pred vs Measured dQ/ds (L_full = split-L)")
    ax_sc.set_xlim(-10, 60)
    ax_sc.set_ylim(-10, 20)
    ax_sc.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    ax_lv.set_xlabel("Timestep")
    ax_lv.set_ylabel(r"H/L_v [-]")
    ax_lv.set_title("All Models: H/L_v")
    ax_lv.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    ax_lk.set_xlabel("Timestep")
    ax_lk.set_ylabel(r"H/L_K [-]")
    ax_lk.set_title("All Models: H/L_K")
    ax_lk.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    ax_lt.set_xlabel("Timestep")
    ax_lt.set_ylabel(r"H/L [-]")
    if const_l_km is not None:
        ax_lt.set_title(f"All Models: H/L (L_full = {const_l_km:g} km)")
    elif use_2lv_as_full:
        ax_lt.set_title(r"All Models: 2H/L_v")
    elif use_lv_as_full:
        ax_lt.set_title(r"All Models: H/L (L_full = L_v)")
    else:
        ax_lt.set_title(r"All Models: H(2/L_v + 1/L_K)")
    ax_lt.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)

    plt.tight_layout()

    if len(files) == 1:
        stem = os.path.basename(files[0]).replace(".txt", "")
    else:
        stem = "multi-model"
    if const_l_km is not None:
        stem = f"{stem}.Lfull-const-{const_l_km:g}km"
    elif use_2lv_as_full:
        stem = f"{stem}.Lfull-eq-2Lv"
    elif use_lv_as_full:
        stem = f"{stem}.Lfull-eq-Lv"
    out_png = os.path.join(out_dir, f"dqds_splitL_compare.{stem}.png")
    out_pdf = os.path.join(out_dir, f"dqds_splitL_compare.{stem}.pdf")

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    main()
