# Session Notes (2026-02-13)

## Repo and setup
- Created GitHub repo: `https://github.com/AdamFHolt/slab-dip-vs-pressure`
- Initialized local git and pushed changes.
- Current branch: `master`.

## High-level project summary
- Root folders:
  - `analysis/` for ASPECT post-processing and plotting.
  - `input_geometries/` for geometry generation scripts.
  - `rheology_params/` for rheology prefactor scripts.
- Data is large (analysis outputs are heavy); version control is configured to avoid tracking huge generated outputs.

## Files added
- `README.md`
- `analysis/ANALYSIS_QUICKSTART.md`
- `analysis/plot_dqds_splitL.py`
- `analysis/many_property-extractions.sh`
- `SESSION_NOTES.md` (this file)

## Key renames/cleanup
- `analysis/extract_propertiesB.py` -> `analysis/extract_properties.py`
- Multiple legacy `B/C/D` names in `analysis/other-scripts/` were renamed to clearer names (`variant`, `testc`, `testd`, `tmp`, etc.).
- `analysis/observations/many_extractionsB.sh` -> `analysis/observations/many_extractions_extra.sh`

## Main scientific/code changes
- Updated `analysis/extract_properties.py` to:
  - write outputs to `analysis/text_files/TESTC/`.
  - auto-create output directory if missing.
  - keep measured in-slab shear gradient term (`slab_stress_term`) from force differencing.
  - add split-length-scale analytical term:
    - `dQds_scaling_splitL = eta * H * K * v_s * (2/L_v + 1/L_K)`
  - compute and save:
    - `vs_mid`, `dvs_ds`, `Lv`, `Lk`, `dQds_scaling_splitL`, `slab_visc_mid`, `slab_visc_x_km`, `slab_visc_y_km`, `L_total_inv=(2/L_v+1/L_K)`
- Output table from `extract_properties.py` is now 30 columns.

## New plotting utility
- `analysis/plot_dqds_splitL.py`:
  - input: one or more `text_files/TESTC/*.txt` outputs from `extract_properties.py`.
  - plots measured vs split-L predicted `dQ/ds` as:
    - time series
    - predicted-vs-measured scatter with 1:1 line and correlation.
  - output location:
    - `analysis/plots/DP-comparisons/compilations/`

## Batch runner
- `analysis/many_property-extractions.sh` now runs `extract_properties.py` in parallel.
- Uses concurrency cap:
  - `MAX_JOBS` env var (default `8`).
- Example:
  - `MAX_JOBS=8 ./many_property-extractions.sh`

## Runtime/environment notes
- Confirmed `analysis` path is on NVMe:
  - `/dev/nvme0n1p9` mounted at `/home` (ext4).
- Desktop CPU: 52 logical CPUs.
- Recommended starting parallelism for this workflow: `MAX_JOBS=8` then test `12`.

## Useful run commands
- Single model extraction:
  - `cd analysis`
  - `python3 extract_properties.py <model_name> <max_time> <analysis_depth_m> <analysis_depth_dz_m> <ds_m> <profile_dz_m>`
- Batch extraction:
  - `cd analysis`
  - `MAX_JOBS=8 ./many_property-extractions.sh`
- Plot split-L comparison:
  - `cd analysis`
  - `python3 plot_dqds_splitL.py text_files/TESTC/<output_file>.txt`

## Recent commits pushed
- `48b3608` refactor: normalize legacy analysis script names
- `82cdd3f` feat: add split-length dQ/ds diagnostics and parallel extraction runner

## Current status
- Work is pushed to GitHub.
- Ready for running new TESTC extractions and evaluating split-L `dQ/ds` scaling.
