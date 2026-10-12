#!/bin/bash
# Extract the variable-plate-thickness runs (H = 70, 130 km x eta' = 250, 500,
# 1000; free, fixed-SP and fixed-OP plates, plus the B-matched drho71 run) at
# one depth, into text_files/extracted_varH/ and without the per-timestep
# evolution figures, so nothing on disk is overwritten.  Models whose output
# file already exists are skipped unless FORCE=1, so new runs can be appended
# to CASES and the script rerun as they come in.
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
#   FORCE=1 bash many_property-extractions.varH.sh     # redo existing outputs too

cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis

OUTDIR="extracted_varH"
SMOOTH="legacy"
DEPTH="${1:-300.0e3}"
MAX_JOBS="${MAX_JOBS:-6}"

ZKM=$(python3 -c "print('%.1f' % (float('$DEPTH')/1e3))")
LOGDIR="text_files/${OUTDIR}/logs"
mkdir -p "$LOGDIR"
echo "outdir = text_files/${OUTDIR}, smoothing = ${SMOOTH}, depth = ${DEPTH}"

# model maxt [drho]   (last output: time, tip depth); drho defaults to 50
CASES=(
"2D_compositional_subd_lower-res_new_250plates_H70km 39"     # 38: 74.7 Myr, 1256 km
"2D_compositional_subd_lower-res_new_250plates_H130km 26"    # 25: 39.9 Myr, 1270 km
"2D_compositional_subd_lower-res_new_500plates_H70km 32"     # 31: 70.0 Myr, 1255 km
"2D_compositional_subd_lower-res_new_500plates_H130km 30"    # 29: 48.7 Myr, 1269 km
"2D_compositional_subd_lower-res_new_1000plates_H70km 31"    # 30: 72.8 Myr, 1262 km
"2D_compositional_subd_lower-res_new_1000plates_H130km 33"   # 32: 53.9 Myr, 1253 km
# fixed OP (added 2026-10-08; runs reached 90 Myr)
"2D_compositional_subd_lower-res_new_FixedOP_250plates_H70km 36"   # 35: 67.5 Myr, 1255 km
"2D_compositional_subd_lower-res_new_FixedOP_500plates_H70km 33"   # 32: 69.7 Myr, 1257 km
"2D_compositional_subd_lower-res_new_FixedOP_1000plates_H70km 36"  # 35: 77.3 Myr, 1250 km
# fixed SP: the tip never reaches 1250 km by 90 Myr, so the window runs to the
# last output (the H = 100 fixed-SP analyses end at 1195-1291 km).  Their
# trenches retreat past x = 2500 km, the left edge of the script's old fixed
# high-res window, which clipped the midplane contour late in these runs
# (shallowest contour depth climbing above 110 km, or no midplane at all for
# the H130 pair); extract_properties.py widens the window itself since
# 2026-10-09 and the four affected models here were re-extracted with FORCE=1.
"2D_compositional_subd_lower-res_new_FixedSP_250plates_H70km 34"   # 33: 88.2 Myr, 1168 km
# added 2026-10-09 (all reached 90 Myr)
"2D_compositional_subd_lower-res_new_FixedOP_250plates_H130km 33"   # 32: 47.8 Myr, 1255 km
"2D_compositional_subd_lower-res_new_FixedOP_500plates_H130km 40"   # 39: 57.4 Myr, 1257 km
"2D_compositional_subd_lower-res_new_FixedOP_1000plates_H130km 42"  # 41: 62.7 Myr, 1253 km
"2D_compositional_subd_lower-res_new_500plates_H70km_drho71 39 71.43"  # 38: 59.8 Myr, 1250 km; B matched to H = 100
# fixed SP, tip short of 1250 km at 90 Myr: window runs to the last output
"2D_compositional_subd_lower-res_new_FixedSP_250plates_H130km 46"   # 45: 89.0 Myr, 1243 km
"2D_compositional_subd_lower-res_new_FixedSP_500plates_H130km 42"   # 41: 87.3 Myr, 1177 km
"2D_compositional_subd_lower-res_new_FixedSP_500plates_H70km 32"    # 31: 86.9 Myr, 1063 km
"2D_compositional_subd_lower-res_new_FixedSP_1000plates_H70km 30"   # 29: 87.5 Myr,  977 km
# resubmitted job 3580803, hit the 48 h limit at 82.4 Myr (synced 2026-10-12)
"2D_compositional_subd_lower-res_new_FixedSP_1000plates_H130km 37"  # 36: 81.1 Myr, 1141.3 km
)

for c in "${CASES[@]}"; do
  set -- $c
  m="$1"; maxt="$2"; drho="${3:-50}"
  if [ -z "${FORCE:-}" ] && [ -e "text_files/${OUTDIR}/$m.z$ZKM.shear-dz10.0.ds10.0.prof-dz1.0km.txt" ]; then
    echo "skipping $m (already extracted at z$ZKM)"; continue
  fi
  PYTHONNOUSERSITE=1 python3 -W ignore extract_properties.py \
      "$m" "$maxt" "$DEPTH" 10.0e3 10.0e3 1.0e3 "$OUTDIR" "$SMOOTH" noplots "$drho" \
      > "$LOGDIR/$m.z$ZKM.log" 2>&1 &
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done
done
wait
echo "ALL EXTRACTIONS DONE (${OUTDIR}, z${ZKM})"
grep -h "^WARNING" "$LOGDIR"/*.z$ZKM.log | wc -l | xargs echo "timesteps flagged as not covered:"
