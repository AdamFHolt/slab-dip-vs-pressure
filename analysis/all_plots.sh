#!/bin/bash
# All plots (from Fig. 3 onward)
# Run from: /home/holt/Projects/ASPECT/subd_2D/compositional/analysis/
# Standard params: analysis_depth=300km, shear-dz=10km, ds=10km, prof-dz=1km

# --- Main figures ---

python3 plot_DPvsDP.color-points.no-ot.py 10000 10000 1000 other

python3 plot_DPvsDP.no-ot.py 10000 10000 1000 other

python3 plot_forces.no-ot.py 300000 10000 10000 1000

python3 plot_forces-bars.no-norm.py 300000 10000 10000 1000

# --- Supplementary figures ---

python3 plot_forces.supp-shear.no-ot.py 300000 10000 10000 1000

python3 plot_Leff_distribution.supp.py

python3 plot_Psp-DP-ratio.py 300000 10000 10000 1000

# --- Single-timestep snapshots (run manually with specific model/timestep) ---
# python3 plot_onestep_simple-pressure.zoomed.py <model_name> <timestep> <x_center_km>
# python3 plot_onestep_simple-viscosity.zoomed.py <model_name> <timestep> <x_center_km>

