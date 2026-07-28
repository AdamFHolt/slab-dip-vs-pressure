#!/bin/bash
# Re-extract the 15-model suite at one depth (300 km by default), into a
# separate output directory and without the per-timestep evolution figures, so
# that nothing already on disk is overwritten.
#
# This exists to test the midplane-coverage fix described in
# extract_properties.py: when a slab necks, the c_llith = 0.5 contour stops
# existing in the shallow section, and dip and curvature "at 300 km" silently
# become edge values read off a contour that starts far deeper.  The new
# columns 30-38 record that, and the smoothing window can now be set as an arc
# length rather than as a fixed number of contour points.
#
# Usage:
#   bash many_property-extractions.recheck.sh TESTD_legacy legacy
#   bash many_property-extractions.recheck.sh TESTD_legacy legacy 250.0e3
#   bash many_property-extractions.recheck.sh TESTD_arc150 150
#
# The original many_property-extractions.sh is unchanged and still writes
# text_files/TESTC/ with the figures.

cd /home/holt/Projects/ASPECT/subd_2D/compositional/analysis

OUTDIR="${1:-TESTD_legacy}"
SMOOTH="${2:-legacy}"
DEPTH="${3:-300.0e3}"
MAX_JOBS="${MAX_JOBS:-15}"

# depth tag for the log names, so several depths can share one output
# directory without overwriting each other's logs (the data files already
# carry the depth in their own filenames)
ZKM=$(python3 -c "print('%.1f' % (float('$DEPTH')/1e3))")
LOGDIR="text_files/${OUTDIR}/logs"
mkdir -p "$LOGDIR"
echo "outdir = text_files/${OUTDIR}, smoothing = ${SMOOTH}, depth = ${DEPTH}"

# model maxt, exactly as in many_property-extractions.sh
CASES=(
"2D_compositional_subd_lower-res_new_50plates 25"
"2D_compositional_subd_lower-res_new_FixedSP_50plates 24"
"2D_compositional_subd_lower-res_new_FixedOP_50plates 25"
"2D_compositional_subd_lower-res_new_250plates 28"
"2D_compositional_subd_lower-res_new_FixedSP_250plates 30"
"2D_compositional_subd_lower-res_new_FixedOP_250plates 28"
"2D_compositional_subd_lower-res_new2 31"
"2D_compositional_subd_FixedSP_lower-res_new2 30"
"2D_compositional_subd_FixedOP_lower-res_new 28"
"2D_compositional_subd_lower-res_new_1000plates 30"
"2D_compositional_subd_lower-res_new_FixedSP_1000plates2 30"
"2D_compositional_subd_lower-res_new_FixedOP_1000plates 30"
"2D_compositional_subd_lower-res_new_375plates 30"
"2D_compositional_subd_lower-res_new_FixedSP_375plates 27"
"2D_compositional_subd_lower-res_new_FixedOP_375plates 32"
)

# SKIP_MODELS is a space-separated list, for when one model is already running
for c in "${CASES[@]}"; do
  set -- $c
  m="$1"; maxt="$2"
  case " ${SKIP_MODELS:-} " in *" $m "*) echo "skipping $m"; continue;; esac
  PYTHONNOUSERSITE=1 python3 -W ignore extract_properties.py \
      "$m" "$maxt" "$DEPTH" 10.0e3 10.0e3 1.0e3 "$OUTDIR" "$SMOOTH" noplots \
      > "$LOGDIR/$m.z$ZKM.log" 2>&1 &
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do wait -n; done
done
wait
echo "ALL EXTRACTIONS DONE (${OUTDIR}, z${ZKM})"
grep -h "^WARNING" "$LOGDIR"/*.z$ZKM.log | wc -l | xargs echo "timesteps flagged as not covered:"
