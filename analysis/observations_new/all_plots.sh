#!/bin/bash

# Log full run output (incl. print_stats summary) to all_plots_run.log (gitignored)
exec > >(tee all_plots_run.log) 2>&1

# H, curvature, vc maps
python3 plot_H-curvature-vc_maps.py

# dip, age maps
python3 plot_dip-age_maps.py

# Lambda breakdown: stress scaling (numerator) and B (denominator) maps
python3 plot_scaling-B_maps.py

# parameter exploration
python3 plot_param_exploration.py

# full map: ref
python3 plot_final_map.py 4e22 3.28e-5 1333 8.044e-7 88e3

# just-maps: ref, and the high-viscosity supplement case (2e23: the first
# viscosity at which a substantial share of segments fail Lambda < 0.1 under
# the cubic scaling; the 2e22/8e22 pair was dropped 2026-10-19)
python3 plot_final_just-maps.py 4e22 3.28e-5 1333 8.044e-7 88e3
python3 plot_final_just-maps.py 2e23 3.28e-5 1333 8.044e-7 88e3

# Lambda / DP summary stats (abstract + Section 3.4 numbers)
python3 print_stats.py

# --- Zip all plot pdfs (fresh zip so removed figures don't linger) ---
rm -f plots/all_observational_plots.zip
zip -j plots/all_observational_plots.zip plots/*.pdf
