# compositional

Tools and data products for 2D compositional subduction experiments (ASPECT).

## Project layout

- `input_geometries/`: Python generators for initial composition and temperature fields used by ASPECT input models.
- `rheology_params/`: Scripts to compute rheological prefactors (diffusion/dislocation style parameter sweeps).
- `analysis/`: Post-processing workflows for ASPECT outputs, including CSV extraction (`pvpython`) and force/stress/dip/curvature diagnostics.

## Typical workflow

1. Generate geometry text files from scripts in `input_geometries/`.
2. Run ASPECT externally with corresponding `.prm` models.
3. Place model outputs under `analysis/raw_outputs/<model_name>/`.
4. Extract per-timestep CSV files with:
   - `pvpython analysis/extract_csv.py <model_name> <max_time>`
5. Compute diagnostics (e.g., `analysis/extract_pref.py`) and plot scripts in `analysis/`.

## Notes

- This repository intentionally ignores very large generated outputs (`analysis/raw_outputs`, `analysis/csv_outputs`, plots, and text outputs) so version control tracks the reproducible scripts/configuration only.
- Some scripts assume local ParaView (`pvpython`) and scientific Python stack (`numpy`, `scipy`, `matplotlib`).
