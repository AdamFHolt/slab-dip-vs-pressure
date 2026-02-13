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
python3 extract_properties.py <model_name> <max_time> <analysis_depth_m> <analysis_depth_dz_m> <ds_m> <profile_dz_m>
```

Example:

```bash
python3 extract_properties.py 2D_compositional_subd_lower-res_new_375plates 33 300000 100000 50000 50000
```

Output:
- `text_files/...` (diagnostic text products)

4. (Optional) Extract far-field pressure time series:

```bash
python3 extract_pref.py <model_name> <max_time> <analysis_depth_m> <analysis_depth_dz_m> <ds_m> <profile_dz_m>
```

5. Generate plots (pick scripts as needed):

```bash
python3 plot_DPvsDP.no-ot.py
python3 plot_forces.no-ot.py
python3 plot_onestep_simple-pressure.zoomed.py <model_name> <timestep>
python3 plot_onestep_simple-viscosity.zoomed.py <model_name> <timestep>
```

Outputs:
- `plots/...`

## Key Script Roles

- `extract_csv.py`: ParaView export from `.pvd` to CSV.
- `extract_properties.py`: main property extraction pipeline.
- `extract_pref.py`: pressure reference extraction.
- `functions.py`: geometry/stress/dip/curvature utility functions.
- `functions_plotting.py`: shared plotting helpers.

## Notes

- Many plotting/extraction scripts assume existing folder names used in this repo.
- If you add a new model name, keep it consistent across `raw_outputs/`, `csv_outputs/`, `text_files/`, and plot commands.
