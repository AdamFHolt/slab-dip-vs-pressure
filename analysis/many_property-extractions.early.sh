#!/bin/bash
# Reviewer 2, comment 3: extract the early outputs (t = 5-10, slab tip ~350-680 km)
# of every model in Figs 3-5 at one depth, into text_files/extracted_early/.
# t = 8-10 are re-extracted too, with the current (window-fixed) code: the H = 130
# extractions of t = 8-10 in extracted_varH predate the fix.
# Same arguments as many_property-extractions.varH.sh, plus first_time = 5.
#
# Usage:  bash many_property-extractions.early.sh [depth_m]   (default 300.0e3)
cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis

OUTDIR="extracted_early"
DEPTH="${1:-300.0e3}"
MAX_JOBS="${MAX_JOBS:-12}"
ZKM=$(python3 -c "print('%.1f' % (float('$DEPTH')/1e3))")
LOGDIR="text_files/${OUTDIR}/logs"
mkdir -p "$LOGDIR"

MODELS=(
2D_compositional_subd_lower-res_new_50plates
2D_compositional_subd_lower-res_new_250plates
2D_compositional_subd_lower-res_new_375plates
2D_compositional_subd_lower-res_new2
2D_compositional_subd_lower-res_new_1000plates
2D_compositional_subd_lower-res_new_FixedSP_50plates
2D_compositional_subd_lower-res_new_FixedSP_250plates
2D_compositional_subd_lower-res_new_FixedSP_375plates
2D_compositional_subd_FixedSP_lower-res_new2
2D_compositional_subd_lower-res_new_FixedSP_1000plates2
2D_compositional_subd_lower-res_new_FixedOP_50plates
2D_compositional_subd_lower-res_new_FixedOP_250plates
2D_compositional_subd_lower-res_new_FixedOP_375plates
2D_compositional_subd_FixedOP_lower-res_new
2D_compositional_subd_lower-res_new_FixedOP_1000plates
)
for bc in "" FixedSP_ FixedOP_; do
  for eta in 250 500 1000; do
    for H in 70 130; do
      MODELS+=("2D_compositional_subd_lower-res_new_${bc}${eta}plates_H${H}km")
    done
  done
done

for m in "${MODELS[@]}"; do
  PYTHONNOUSERSITE=1 python3 -W ignore extract_properties.py \
      "$m" 11 "$DEPTH" 10.0e3 10.0e3 1.0e3 "$OUTDIR" legacy noplots 50 5 \
      > "$LOGDIR/$m.z$ZKM.log" 2>&1 &
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done
done
wait
echo "ALL EARLY EXTRACTIONS DONE (${OUTDIR}, z${ZKM}): $(ls text_files/${OUTDIR}/*.z$ZKM.*.txt | wc -l) files"
