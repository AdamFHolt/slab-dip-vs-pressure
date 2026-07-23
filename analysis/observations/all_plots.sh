#!/bin/bash

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

# just-maps: ref, low visc, high visc
python3 plot_final_just-maps.py 4e22 3.28e-5 1333 8.044e-7 88e3
python3 plot_final_just-maps.py 2e22 3.28e-5 1333 8.044e-7 88e3
python3 plot_final_just-maps.py 8e22 3.28e-5 1333 8.044e-7 88e3

# Lambda / DP summary stats (abstract + Section 3.4 numbers)
python3 print_stats.py

# --- Zip all plot pdfs (fresh zip so removed figures don't linger) ---
rm -f plots/all_observational_plots.zip
zip -j plots/all_observational_plots.zip plots/*.pdf
