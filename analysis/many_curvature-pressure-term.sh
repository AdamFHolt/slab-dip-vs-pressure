#!/bin/bash
# Run compute_curvature_pressure_term.py over the 15 manuscript models (300 km).
cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis
MODELS=(
2D_compositional_subd_lower-res_new_50plates
2D_compositional_subd_lower-res_new_FixedSP_50plates
2D_compositional_subd_lower-res_new_FixedOP_50plates
2D_compositional_subd_lower-res_new_250plates
2D_compositional_subd_lower-res_new_FixedSP_250plates
2D_compositional_subd_lower-res_new_FixedOP_250plates
2D_compositional_subd_lower-res_new_375plates
2D_compositional_subd_lower-res_new_FixedSP_375plates
2D_compositional_subd_lower-res_new_FixedOP_375plates
2D_compositional_subd_lower-res_new2
2D_compositional_subd_FixedSP_lower-res_new2
2D_compositional_subd_FixedOP_lower-res_new
2D_compositional_subd_lower-res_new_1000plates
2D_compositional_subd_lower-res_new_FixedSP_1000plates2
2D_compositional_subd_lower-res_new_FixedOP_1000plates
)
mkdir -p text_files/curvature_pressure_term/logs
for m in "${MODELS[@]}"; do
  PYTHONNOUSERSITE=1 python3 -W ignore compute_curvature_pressure_term.py "$m" \
    > text_files/curvature_pressure_term/logs/"$m".log 2>&1 &
done
wait
echo "ALL DONE"
