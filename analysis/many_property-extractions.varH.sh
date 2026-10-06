#!/bin/bash
# Extract the six variable-plate-thickness runs (H = 70, 130 km x eta' = 250,
# 500, 1000; free plates) at one depth, into text_files/extracted_varH/ and
# without the per-timestep evolution figures, so nothing on disk is overwritten.
#
# Arguments match many_property-extractions.recheck.sh, whose extracted_coverage/
# output holds the matching H = 100 runs (new_250plates, new2, new_1000plates)
# from the same version of extract_properties.py.
#
# maxt: the first output at which the slab tip (deepest c_llith > 0.5) passes
# 1250 km is the last one analysed (range() is exclusive, hence +1).  The
# H = 100 analyses end with the tip at 1195-1291 km.  Output 11, the first one
# the figures use, falls at 642-719 km in all six runs (646-709 km at H = 100).
#
# Usage:
#   bash many_property-extractions.varH.sh               # 300 km
#   bash many_property-extractions.varH.sh 250.0e3

cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis

OUTDIR="extracted_varH"
SMOOTH="legacy"
DEPTH="${1:-300.0e3}"
MAX_JOBS="${MAX_JOBS:-6}"

ZKM=$(python3 -c "print('%.1f' % (float('$DEPTH')/1e3))")
LOGDIR="text_files/${OUTDIR}/logs"
mkdir -p "$LOGDIR"
echo "outdir = text_files/${OUTDIR}, smoothing = ${SMOOTH}, depth = ${DEPTH}"

# model maxt   (last output: time, tip depth)
CASES=(
"2D_compositional_subd_lower-res_new_250plates_H70km 39"     # 38: 74.7 Myr, 1256 km
"2D_compositional_subd_lower-res_new_250plates_H130km 26"    # 25: 39.9 Myr, 1270 km
"2D_compositional_subd_lower-res_new_500plates_H70km 32"     # 31: 70.0 Myr, 1255 km
"2D_compositional_subd_lower-res_new_500plates_H130km 30"    # 29: 48.7 Myr, 1269 km
"2D_compositional_subd_lower-res_new_1000plates_H70km 31"    # 30: 72.8 Myr, 1262 km
"2D_compositional_subd_lower-res_new_1000plates_H130km 33"   # 32: 53.9 Myr, 1253 km
)

for c in "${CASES[@]}"; do
  set -- $c
  m="$1"; maxt="$2"
  PYTHONNOUSERSITE=1 python3 -W ignore extract_properties.py \
      "$m" "$maxt" "$DEPTH" 10.0e3 10.0e3 1.0e3 "$OUTDIR" "$SMOOTH" noplots \
      > "$LOGDIR/$m.z$ZKM.log" 2>&1 &
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done
done
wait
echo "ALL EXTRACTIONS DONE (${OUTDIR}, z${ZKM})"
grep -h "^WARNING" "$LOGDIR"/*.z$ZKM.log | wc -l | xargs echo "timesteps flagged as not covered:"
