# Analysis Quickstart

This folder processes ASPECT model output into CSVs, derived metrics, and plots.

## Prerequisites

- Python with: `numpy`, `scipy`, `matplotlib`
- ParaView `pvpython` (for `extract_csv.py`)
- Model output directory:
  - `raw_outputs/<model_name>/solution.pvd`
  - `raw_outputs/<model_name>/statistics`

## Recommended Run Order

1. Go to this folder:

```bash
cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis
```

2. Extract per-timestep CSV files from ParaView output:

```bash
pvpython extract_csv.py <model_name> <max_time>
```

Example:

```bash
pvpython extract_csv.py 2D_compositional_subd_lower-res_new_375plates 33
```

Output:
- `csv_outputs/<model_name>/full.<timestep>.csv`

3. Extract diagnostics used by later plots:

```bash
python3 extract_properties.py <model_name> <max_time> <analysis_depth_m> <analysis_depth_dz_m> <ds_m> <profile_dz_m> [outdir] [smoothing] [noplots]
```

Example, using the manuscript's standard parameters:

```bash
python3 extract_properties.py 2D_compositional_subd_lower-res_new_375plates 33 300000 10000 10000 1000
```

Output:
- `text_files/extracted_Lscales/<model_name>.z<depth>...txt`

Whole suite at once: `bash many_property-extractions.sh`.

The three optional arguments default to the historical behaviour: output
subdirectory (default `extracted_Lscales`), Savitzky-Golay smoothing
(`legacy` for 601 contour points, or a length in km), and `noplots` to skip
the per-timestep evolution figures, which dominate the runtime.

4. (Optional) Extract far-field pressure time series:

```bash
python3 extract_pref.py <model_name> <max_time> <analysis_depth_m> <analysis_depth_dz_m> <ds_m> <profile_dz_m>
```

5. Measure the curvature-pressure term and fold it into the normal-stress
column. Every force figure reads the result, so this must be run before them:

```bash
for z in 250.0e3 300.0e3 350.0e3; do bash many_curvature-pressure-term.sh $z; done
python3 make_withT.py
```

6. (Optional) Flag timesteps where the slab midplane does not resolve the
analysis depth, and blank them:

```bash
for z in 250.0e3 300.0e3 350.0e3; do
  bash many_property-extractions.recheck.sh extracted_coverage legacy $z
done
python3 make_withT_covered.py
```

7. Generate the manuscript figures:

```bash
bash all_plots.sh
```

That runs `make_withT.py`, the seven figure scripts, and bundles the eight
PDFs into `plots/DP-comparisons/compilations/all_model_plots.zip`. Individual
scripts take the same parameters, for example:

```bash
python3 plot_DPvsDP.no-ot.py 10000 10000 1000 other
python3 plot_forces.no-ot.py 300000 10000 10000 1000
python3 plot_onestep_simple-pressure.zoomed.py <model_name> <timestep> <x_center_km>
```

The coverage-filtered variants are `plot_DPvsDP.no-ot.covered.py` and
`plot_DPvsDP.color-points.no-ot.covered.py`, same arguments.

Outputs:
- `plots/...`

## Key Script Roles

- `extract_csv.py`: ParaView export from `.pvd` to CSV.
- `extract_properties.py`: main property extraction pipeline.
- `extract_pref.py`: pressure reference extraction.
- `functions.py`: geometry/stress/dip/curvature utility functions.
- `functions_plotting.py`: shared plotting helpers.
- `compute_curvature_pressure_term.py`: measures the pressure part of the
  curvature term, `T = (H/2) K (2*Pbar_n - P_subslab - P_wedge)`.
- `make_withT.py`: adds `T` to column 17 of the archived extraction.
- `make_withT_covered.py`: blanks timesteps failing the midplane coverage test.
- `all_plots.sh`: the manuscript figures, end to end.

`text_files/README.txt` describes what each data directory holds and which
script writes and reads it.

## Notes

- Timesteps start at `t = 8` in the extraction; the figures plot from `t = 11`
  (`tactual_min = 11`, `tmin = 3`), because the earliest timesteps close the
  force balance far less well while the slab is still establishing itself.
- Many plotting/extraction scripts assume existing folder names used in this repo.
- If you add a new model name, keep it consistent across `raw_outputs/`,
  `csv_outputs/`, `text_files/`, and plot commands. The model list is repeated
  in `make_withT.py`, `make_withT_covered.py`, and the `many_*.sh` runners.
