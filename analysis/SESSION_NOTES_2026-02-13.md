# Session Notes (2026-02-13)

## Summary
- Updated `plot_dqds_splitL.py` from a 2-panel plot to a 4-panel diagnostic figure.
- Removed legends from all panels.
- Switched panel 1 scatter to:
  - x: `Predicted dQ/ds [MPa]`
  - y: `-Measured dQ/ds [MPa]`
- Added zoom to panel 1:
  - x-limits: `[-10, 60]` MPa
  - y-limits: `[-10, 20]` MPa
- Changed lower panels to non-dimensional forms:
  - `H/L_v`
  - `H/L_K`
  - `H/L` (depends on chosen mode)

## Plot Modes Added
- Default mode (no flag):
  - Uses saved split-length term from file: `1/L = 2/L_v + 1/L_K`
- `--full-L-equals-Lv`
  - Uses `1/L = 1/L_v`
- `--full-L-equals-2Lv`
  - Uses `1/L = 2/L_v`
- `--const-L-km <value>`
  - Uses constant full length: `1/L = 1/(value_km * 1000)`

These modes are mutually exclusive.

## Current Script Behavior
- If no input files are provided, script auto-loads:
  - `text_files/TESTC/*.txt`
- Prediction is recomputed in-plot from components using:
  - `pred_pa = eta * H * (K * (vs / seconds_per_year)) * L_inv`
  - converted to MPa for plotting.

## Files Generated During Session
- `plots/DP-comparisons/compilations/dqds_splitL_compare.multi-model.Lfull-const-1000km.png`
- `plots/DP-comparisons/compilations/dqds_splitL_compare.multi-model.Lfull-const-1000km.pdf`
- `plots/DP-comparisons/compilations/dqds_splitL_compare.multi-model.Lfull-eq-Lv.png`
- `plots/DP-comparisons/compilations/dqds_splitL_compare.multi-model.Lfull-eq-Lv.pdf`
- `plots/DP-comparisons/compilations/dqds_splitL_compare.multi-model.Lfull-eq-2Lv.png`
- `plots/DP-comparisons/compilations/dqds_splitL_compare.multi-model.Lfull-eq-2Lv.pdf`
- (plus other intermediate variants in same folder)

## Diagnostic Findings
- `1/L_K` strongly dominates `2/L_v + 1/L_K` in most timesteps/models.
- Quick aggregate check showed:
  - median contribution from `1/L_K` to total `L_inv` ~ `0.98`
  - majority of timesteps had `> 80%` contribution from `1/L_K`
- Interpretation: overprediction is likely tied to `L_K = |K / (dK/ds)|` behavior and/or `dK/ds` sensitivity.

## Suggested Next Steps
1. Inspect/fix `L_K` pipeline in `functions.py` / `extract_properties.py`.
2. Add safeguards for unstable `K/(dK/ds)` (e.g., `|K|` floor, robust derivative).
3. Re-run extraction for a small subset and compare old vs new `L_K`, `L_inv`, and predicted dQ/ds.

## Priority For Next Session
- Write a highly pared-down extractor focused only on length scales:
  - proposed file: `extract_L_scales.py`
  - goal: carefully compute and save only `L_v` and `L_K` (+ minimal raw audit terms)
  - include explicit units and minimal dependencies
  - start with one/few models for verification before broader reruns
