#!/bin/bash
# Curvature-pressure term T for the variable-H comparison set: the six
# H = 70/130 km runs and their three H = 100 counterparts.  Writes
# text_files/curvature_pressure_term_varH/ only; see
# compute_curvature_pressure_term.varH.py.
# Usage: bash many_curvature-pressure-term.varH.sh [analysis_depth_m]   (default 300.0e3)
cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis
DEPTH="${1:-300.0e3}"
ZKM=$(python3 -c "print('%.1f' % (float('$DEPTH')/1e3))")
echo "analysis depth = $ZKM km"
MODELS=(
2D_compositional_subd_lower-res_new_250plates_H70km
2D_compositional_subd_lower-res_new_250plates
2D_compositional_subd_lower-res_new_250plates_H130km
2D_compositional_subd_lower-res_new_500plates_H70km
2D_compositional_subd_lower-res_new2
2D_compositional_subd_lower-res_new_500plates_H130km
2D_compositional_subd_lower-res_new_1000plates_H70km
2D_compositional_subd_lower-res_new_1000plates
2D_compositional_subd_lower-res_new_1000plates_H130km
)
mkdir -p text_files/curvature_pressure_term_varH/logs
for m in "${MODELS[@]}"; do
  PYTHONNOUSERSITE=1 python3 -W ignore compute_curvature_pressure_term.varH.py "$m" "$DEPTH" \
    > text_files/curvature_pressure_term_varH/logs/"$m".z"$ZKM".log 2>&1 &
done
wait
echo "ALL DONE"
