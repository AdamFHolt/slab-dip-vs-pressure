#!/bin/bash
# Curvature-pressure term T for the variable-H comparison set: the H = 70/130 km
# runs (free, fixed-SP, fixed-OP, drho71) and their H = 100 counterparts.
# Models whose extraction file is not there yet are skipped, so the list can
# hold the whole set while the runs come in.  Writes
# text_files/curvature_pressure_term_varH/ only; see
# compute_curvature_pressure_term.varH.py.
# Usage: bash many_curvature-pressure-term.varH.sh [analysis_depth_m]   (default 300.0e3)
cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis
DEPTH="${1:-300.0e3}"
ZKM=$(python3 -c "print('%.1f' % (float('$DEPTH')/1e3))")
echo "analysis depth = $ZKM km"
MODELS=(
# free plates
2D_compositional_subd_lower-res_new_250plates_H70km
2D_compositional_subd_lower-res_new_250plates
2D_compositional_subd_lower-res_new_250plates_H130km
2D_compositional_subd_lower-res_new_500plates_H70km
2D_compositional_subd_lower-res_new2
2D_compositional_subd_lower-res_new_500plates_H130km
2D_compositional_subd_lower-res_new_1000plates_H70km
2D_compositional_subd_lower-res_new_1000plates
2D_compositional_subd_lower-res_new_1000plates_H130km
# free, B matched to H = 100 (drho = 71.43)
2D_compositional_subd_lower-res_new_500plates_H70km_drho71
# fixed SP
2D_compositional_subd_lower-res_new_FixedSP_250plates_H70km
2D_compositional_subd_lower-res_new_FixedSP_250plates
2D_compositional_subd_lower-res_new_FixedSP_250plates_H130km
2D_compositional_subd_lower-res_new_FixedSP_500plates_H70km
2D_compositional_subd_FixedSP_lower-res_new2
2D_compositional_subd_lower-res_new_FixedSP_500plates_H130km
2D_compositional_subd_lower-res_new_FixedSP_1000plates_H70km
2D_compositional_subd_lower-res_new_FixedSP_1000plates2
2D_compositional_subd_lower-res_new_FixedSP_1000plates_H130km
# fixed OP
2D_compositional_subd_lower-res_new_FixedOP_250plates_H70km
2D_compositional_subd_lower-res_new_FixedOP_250plates
2D_compositional_subd_lower-res_new_FixedOP_250plates_H130km
2D_compositional_subd_lower-res_new_FixedOP_500plates_H70km
2D_compositional_subd_FixedOP_lower-res_new
2D_compositional_subd_lower-res_new_FixedOP_500plates_H130km
2D_compositional_subd_lower-res_new_FixedOP_1000plates_H70km
2D_compositional_subd_lower-res_new_FixedOP_1000plates
2D_compositional_subd_lower-res_new_FixedOP_1000plates_H130km
)
SUFFIX="shear-dz10.0.ds10.0.prof-dz1.0km.txt"
mkdir -p text_files/curvature_pressure_term_varH/logs
for m in "${MODELS[@]}"; do
  if [ ! -e text_files/extracted_varH/"$m".z"$ZKM"."$SUFFIX" ] && [ ! -e text_files/extracted_coverage/"$m".z"$ZKM"."$SUFFIX" ]; then
    echo "skipping $m (no extraction at z$ZKM yet)"; continue
  fi
  PYTHONNOUSERSITE=1 python3 -W ignore compute_curvature_pressure_term.varH.py "$m" "$DEPTH" \
    > text_files/curvature_pressure_term_varH/logs/"$m".z"$ZKM".log 2>&1 &
done
wait
echo "ALL DONE"
